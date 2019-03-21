
IF(BUILD_LIBADC)
    
    SET( LIBADC1_SOURCES   src/sizes.F KDTREE2/kdtree2.F 
                           src/global.F src/boundaries.F src/global_3dvs.F
                           src/messenger.F )

    SET( LIBADC2_SOURCES   src/mesh.F src/harm.F wind/vortex.F src/wind.F src/hashtable.F
                           src/owiwind.F src/rs2.F src/owi_ice.F 
                           src/itpackv.F src/nodalattr.F src/globalio.F 
                           src/subdomain.F src/gwce.F src/wetdry.F src/momentum.F
                           src/netcdfio.F src/control.F src/xdmfio.F )

    SET( LIBADC3_SOURCES    src/writer.F )

    SET( LIBADC4_SOURCES   src/write_output.F src/couple2swan.F src/adcirc.F
                           src/weir_boundary.F src/read_input.F src/cstart.F src/hstart.F 
                           src/timestep.F src/vsmy.F src/transport.F src/driver.F )

    ADD_LIBRARY(templib_libadc1 ${LIBADC1_SOURCES})
    ADD_LIBRARY(templib_libadc2 ${LIBADC2_SOURCES})
    ADD_LIBRARY(templib_libadc3 ${LIBADC3_SOURCES})
    
    ADD_library(adc ${LIBADC4_SOURCES})

    addCompilerFlags(templib_libadc1)
    addCompilerFlags(templib_libadc2)
    addCompilerFlags(templib_libadc3)
    addCompilerFlags(adc)

    addMPI(templib_libadc1)
    addMPI(templib_libadc2)
    addMPI(templib_libadc3)
    addMPI(adc)
    
    addLibMkdir(templib_libadc1)
    addLibMkdir(templib_libadc2)
    addLibMkdir(templib_libadc3)
    addLibMkdir(adc)
    addLibVersion(adc)

    SET_TARGET_PROPERTIES( adc PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )

    TARGET_INCLUDE_DIRECTORIES(templib_libadc2  PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/mod/templib_libadc1)
    TARGET_INCLUDE_DIRECTORIES(templib_libadc3  PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/mod/templib_libadc1)
    TARGET_INCLUDE_DIRECTORIES(adc              PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/mod/templib_libadc1)
    TARGET_INCLUDE_DIRECTORIES(templib_libadc3  PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/mod/templib_libadc2)
    TARGET_INCLUDE_DIRECTORIES(adc              PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/mod/templib_libadc2)
    TARGET_INCLUDE_DIRECTORIES(adc              PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/mod/templib_libadc3)

    TARGET_LINK_LIBRARIES(adc templib_libadc1 templib_libadc2 templib_libadc3)

    ADD_DEPENDENCIES(adc              templib_libadc3)
    ADD_DEPENDENCIES(templib_libadc2  templib_libadc1)
    ADD_DEPENDENCIES(templib_libadc3  templib_libadc2)
    ADD_DEPENDENCIES(templib_libadc1  version)
    ADD_DEPENDENCIES(templib_libadc1  mkdir)
    
    INSTALL(TARGETS adc RUNTIME DESTINATION lib ARCHIVE DESTINATION lib)

ENDIF(BUILD_LIBADC)
