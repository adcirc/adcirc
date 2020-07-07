SET( LIBADC_SOURCES   ${CMAKE_SOURCE_DIR}/src/sizes.F 
                      ${CMAKE_SOURCE_DIR}/thirdparty/KDTREE2/kdtree2.F
                      ${CMAKE_SOURCE_DIR}/src/global.F 
                      ${CMAKE_SOURCE_DIR}/src/boundaries.F 
                      ${CMAKE_SOURCE_DIR}/src/global_3dvs.F
                      ${CMAKE_SOURCE_DIR}/src/messenger.F 
                      ${CMAKE_SOURCE_DIR}/src/mesh.F 
                      ${CMAKE_SOURCE_DIR}/src/harm.F 
                      ${CMAKE_SOURCE_DIR}/wind/vortex.F 
                      ${CMAKE_SOURCE_DIR}/src/wind.F 
                      ${CMAKE_SOURCE_DIR}/src/hashtable.F 
                      ${CMAKE_SOURCE_DIR}/src/owiwind.F 
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
                      ${CMAKE_SOURCE_DIR}/src/quadtrature.F 
                      ${CMAKE_SOURCE_DIR}/src/couple2baroclinic3D.F )

IF(BUILD_LIBADCIRC_STATIC)

    ADD_LIBRARY(libadcirc_static STATIC ${LIBADC_SOURCES})
    INSTALL(TARGETS libadcirc_static ARCHIVE DESTINATION lib)
    ADD_DEPENDENCIES(libadcirc_static version mkdir)
    TARGET_LINK_LIBRARIES(libadcirc_static version mkdir)
    SET_TARGET_PROPERTIES( libadcirc_static PROPERTIES OUTPUT_NAME "adcirc" )
    addCompilerFlags(libadcirc_static)
    addMPI(libadcirc_static)
    SET_TARGET_PROPERTIES( libadcirc_static PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )
    
    INSTALL(TARGETS libadcirc_static RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	                             ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
				     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

ENDIF(BUILD_LIBADCIRC_STATIC)

IF(BUILD_LIBADCIRC_SHARED)
    SET( LIBMKDIR2_SOURCES prep/mkdir.c )
    SET( LIBADC_SHARED_SOURCES  ${LIBADC_SOURCES} ${CMAKE_BINARY_DIR}/CMakeFiles/version_cmake.F )

    ADD_LIBRARY(mkdir2 STATIC ${LIBMKDIR2_SOURCES})
    ADD_LIBRARY(libadcirc_shared SHARED ${LIBADC_SHARED_SOURCES})

    addCompilerFlags(libadcirc_shared)
    addMPI(libadcirc_shared)

    ADD_DEPENDENCIES(libadcirc_shared version mkdir2 )
    TARGET_LINK_LIBRARIES( libadcirc_shared mkdir2 )

    SET_TARGET_PROPERTIES( libadcirc_shared PROPERTIES OUTPUT_NAME "adcirc" )
    
    SET_PROPERTY(TARGET libadcirc_shared PROPERTY POSITION_INDEPENDENT_CODE ON)
    SET_PROPERTY(TARGET mkdir2 PROPERTY POSITION_INDEPENDENT_CODE ON)

    SET_TARGET_PROPERTIES( libadcirc_shared PROPERTIES VERSION ${ADCIRC_VERSION_STRING} SOVERSION ${ADCIRC_VERSION_MAJOR} )
    WRITE_BASIC_PACKAGE_VERSION_FILE( libadcircConfigVersion.cmake VERSION ${ADCIRC_VERSION_STRING} COMPATIBILITY SameMajorVersion )

    INSTALL(TARGETS libadcirc_shared RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	                             ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
				     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
    
    INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/libadcircConfigVersion.cmake DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake )			     

    IF(${CMAKE_INSTALL_PREFIX} STREQUAL "/usr/local")
        INSTALL ( CODE
          "EXECUTE_PROCESS (COMMAND \"${CMAKE_COMMAND}\" -E copy_directory \"${CMAKE_BINARY_DIR}/CMakeFiles/mod/libadcirc_shared\" \"${CMAKE_INSTALL_PREFIX}/include/adcirc\")")
    ELSE()
        INSTALL ( CODE
          "EXECUTE_PROCESS (COMMAND \"${CMAKE_COMMAND}\" -E copy_directory \"${CMAKE_BINARY_DIR}/CMakeFiles/mod/libadcirc_shared\" \"${CMAKE_INSTALL_PREFIX}/include\")")
    ENDIF()

ENDIF(BUILD_LIBADCIRC_SHARED)
