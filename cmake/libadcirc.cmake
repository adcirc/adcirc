SET( LIBADC_SOURCES   src/sizes.F KDTREE2/kdtree2.F 
                      src/global.F src/boundaries.F src/global_3dvs.F
                      src/messenger.F src/mesh.F src/harm.F wind/vortex.F 
                      src/wind.F src/hashtable.F src/owiwind.F src/rs2.F src/owi_ice.F 
                      src/itpackv.F src/nodalattr.F src/globalio.F 
                      src/subdomain.F src/gwce.F src/wetdry.F src/momentum.F
                      src/netcdfio.F src/control.F src/xdmfio.F src/writer.F 
                      src/write_output.F src/couple2swan.F src/adcirc.F
                      src/weir_boundary.F src/read_input.F src/cstart.F src/hstart.F 
                      src/timestep.F src/vsmy.F src/transport.F )

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
