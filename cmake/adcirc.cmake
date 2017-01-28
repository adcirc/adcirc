
IF(BUILD_ADCIRC)
    
    #...Note that we need to build ADCIRC in steps to correctly generate
    #   the objects in the correct order of dependencies. CMAKE tries to
    #   generate this order automatically, however, seems to get confused
    #   by heavy use of compiler flags
    
    SET( ADCIRC1_SOURCES  src/sizes.F KDTREE2/kdtree2.F 
                          src/global.F src/boundaries.F src/mesh.F 
                          src/global_3dvs.F src/harm.F wind/vortex.F 
                          src/wind.F src/owiwind.F src/rs2.F 
                          src/owi_ice.F src/itpackv.F src/nodalattr.F 
                          src/globalio.F src/netcdfio.F src/control.F
                          src/xdmfio.F src/hashtable.F )

    SET( ADCIRC2_SOURCES  src/write_output.F 
                          src/couple2swan.F src/adcirc.F src/subdomain.F 
                          src/weir_boundary.F src/read_input.F src/cstart.F 
                          src/hstart.F src/timestep.F src/vsmy.F 
                          src/transport.F src/driver.F )

    ADD_LIBRARY(templib_adcirc1 ${ADCIRC1_SOURCES})
    
    ADD_EXECUTABLE(adcirc ${ADCIRC2_SOURCES})
    
    SET(ADCIRC_COMPILER_FLAGS "${Fortran_LINELENGTH_FLAG} ${Fortran_COMPILER_SPECIFIC_FLAG} ${PRECISION_FLAG} ${ADCIRC_OPTION_FLAGS} ${ADDITIONAL_FLAGS_ADCIRC}")
    
    IF(NETCDF_WORKING AND NOT XDMF_WORKING)
        SET(ADCIRC_COMPILER_FLAGS "${ADCIRC_COMPILER_FLAGS} ${NETCDF_FLAG} ${NETCDF_COMPRESSION_FLAG}")
        SET_TARGET_PROPERTIES(adcirc PROPERTIES LINK_FLAGS ${NETCDF_LINKER_FLAG})
        TARGET_LINK_LIBRARIES(adcirc templib_adcirc1 netcdf netcdff)
    ELSEIF(NETCDF_WORKING AND XDMF_WORKING)
        SET(ADCIRC_COMPILER_FLAGS "${ADCIRC_COMPILER_FLAGS} ${NETCDF_FLAG} ${NETCDF_COMPRESSION_FLAG} ${XDMF_FLAG}")
        SET(ADCIRC_LINKER_FLAGS "${NETCDF_LINKER_FLAG} ${XDMF_LINKER_FLAG}")
        SET_TARGET_PROPERTIES(adcirc PROPERTIES LINK_FLAGS ${ADCIRC_LINKER_FLAGS} )
        TARGET_LINK_LIBRARIES(adcirc templib_adcirc1 netcdf netcdff XdmfCore XdmfUtils Xdmf)
        TARGET_INCLUDE_DIRECTORIES(adcirc PRIVATE ${CMAKE_SOURCE_DIR}/src)
    ELSE(NETCDF_WORKING AND NOT XDMF_WORKING)
        TARGET_LINK_LIBRARIES(adcirc templib_adcirc1)
    ENDIF(NETCDF_WORKING AND NOT XDMF_WORKING)

    TARGET_INCLUDE_DIRECTORIES(adcirc          PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    TARGET_INCLUDE_DIRECTORIES(templib_adcirc1 PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)

    SET_TARGET_PROPERTIES(templib_adcirc1 PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/adcirc_mod)
    SET_TARGET_PROPERTIES(templib_adcirc1 PROPERTIES COMPILE_FLAGS ${ADCIRC_COMPILER_FLAGS})
    SET_TARGET_PROPERTIES(adcirc PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/adcirc_mod)
    SET_TARGET_PROPERTIES(adcirc PROPERTIES COMPILE_FLAGS ${ADCIRC_COMPILER_FLAGS})
    
    TARGET_LINK_LIBRARIES(adcirc  version)
    TARGET_LINK_LIBRARIES(templib_adcirc1 version)
    
    ADD_DEPENDENCIES(adcirc          templib_adcirc1)
    ADD_DEPENDENCIES(templib_adcirc1 version)
    
    INSTALL(TARGETS adcirc RUNTIME DESTINATION bin)

ENDIF(BUILD_ADCIRC)
