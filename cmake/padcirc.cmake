
IF(BUILD_PADCIRC)
    
    #...Note that we need to build PADCIRC in steps to correctly generate
    #   the objects in the correct order of dependencies. CMAKE tries to
    #   generate this order automatically, however, seems to get confused
    #   by heavy use of compiler flags
    
    SET( PADCIRC1_SOURCES  src/sizes.F KDTREE2/kdtree2.F 
                           src/global.F src/boundaries.F src/global_3dvs.F
                           src/messenger.F )

    SET( PADCIRC2_SOURCES  src/mesh.F src/harm.F wind/vortex.F src/wind.F 
                           src/owiwind.F src/rs2.F src/owi_ice.F 
                           src/itpackv.F src/nodalattr.F src/globalio.F 
                           src/netcdfio.F src/control.F src/xdmfio.F )

    SET( PADCIRC3_SOURCES  src/writer.F )

    SET( PADCIRC4_SOURCES  src/write_output.F src/couple2swan.F src/adcirc.F src/subdomain.F 
                           src/weir_boundary.F src/wetdry.F src/gwce.F src/momentum.F 
                           src/read_input.F src/cstart.F src/hstart.F 
                           src/timestep.F src/vsmy.F src/transport.F src/driver.F )

    ADD_LIBRARY(templib_padcirc1 ${PADCIRC1_SOURCES})
    ADD_LIBRARY(templib_padcirc2 ${PADCIRC2_SOURCES})
    ADD_LIBRARY(templib_padcirc3 ${PADCIRC3_SOURCES})
    
    ADD_EXECUTABLE(padcirc ${PADCIRC4_SOURCES})

    TARGET_INCLUDE_DIRECTORIES(templib_padcirc1 PRIVATE prep)
    TARGET_INCLUDE_DIRECTORIES(padcirc          PRIVATE prep)
    
    SET(PADCIRC_COMPILER_FLAGS "${Fortran_LINELENGTH_FLAG} ${Fortran_COMPILER_SPECIFIC_FLAG} ${PRECISION_FLAG} ${ADCIRC_MPI_FLAG} ${ADCIRC_OPTION_FLAGS} ${ADDITIONAL_FLAGS_ADCIRC}")
    
    IF(NETCDF_WORKING AND NOT XDMF_WORKING)
        SET(PADCIRC_COMPILER_FLAGS "${PADCIRC_COMPILER_FLAGS} ${NETCDF_FLAG} ${NETCDF_COMPRESSION_FLAG}")
        SET_TARGET_PROPERTIES(padcirc PROPERTIES LINK_FLAGS ${NETCDF_LINKER_FLAG})
    ELSEIF(NETCDF_WORKING AND XDMF_WORKING)
        SET(PADCIRC_COMPILER_FLAGS "${PADCIRC_COMPILER_FLAGS} ${NETCDF_FLAG} ${NETCDF_COMPRESSION_FLAG} ${XDMF_FLAG}")
        SET(PADCIRC_LINKER_FLAGS "${NETCDF_LINKER_FLAG} ${XDMF_LINKER_FLAG}")
        SET_TARGET_PROPERTIES(padcirc PROPERTIES LINK_FLAGS ${PADCIRC_LINKER_FLAGS} )
        TARGET_INCLUDE_DIRECTORIES(padcirc PRIVATE ${CMAKE_SOURCE_DIR}/src)
    ENDIF(NETCDF_WORKING AND NOT XDMF_WORKING)
    
    SET_TARGET_PROPERTIES(templib_padcirc1 PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/padcirc_mod)
    SET_TARGET_PROPERTIES(templib_padcirc1 PROPERTIES COMPILE_FLAGS ${PADCIRC_COMPILER_FLAGS})

    SET_TARGET_PROPERTIES(templib_padcirc2 PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/padcirc_mod)
    SET_TARGET_PROPERTIES(templib_padcirc2 PROPERTIES COMPILE_FLAGS ${PADCIRC_COMPILER_FLAGS})
    
    SET_TARGET_PROPERTIES(templib_padcirc3 PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/padcirc_mod)
    SET_TARGET_PROPERTIES(templib_padcirc3 PROPERTIES COMPILE_FLAGS ${PADCIRC_COMPILER_FLAGS})

    SET_TARGET_PROPERTIES(padcirc PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/padcirc_mod)
    SET_TARGET_PROPERTIES(padcirc PROPERTIES COMPILE_FLAGS ${PADCIRC_COMPILER_FLAGS})
    
    TARGET_INCLUDE_DIRECTORIES(templib_padcirc1  PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    TARGET_INCLUDE_DIRECTORIES(templib_padcirc2  PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    TARGET_INCLUDE_DIRECTORIES(templib_padcirc3  PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    TARGET_INCLUDE_DIRECTORIES(padcirc           PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/version_mod)
    
    TARGET_INCLUDE_DIRECTORIES(templib_padcirc1  PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    TARGET_INCLUDE_DIRECTORIES(templib_padcirc2  PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    TARGET_INCLUDE_DIRECTORIES(templib_padcirc3  PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    TARGET_INCLUDE_DIRECTORIES(padcirc           PRIVATE ${MPI_Fortran_INCLUDE_PATH})

    IF(NETCDF_WORKING AND NOT XDMF_WORKING)
        TARGET_LINK_LIBRARIES(padcirc templib_padcirc3 templib_padcirc2 templib_padcirc1 
                                      version mkdir netcdf netcdff ${MPI_Fortran_LIBRARIES})
    ELSEIF(NETCDF_WORKING AND XDMF_WORKING)
        TARGET_LINK_LIBRARIES(padcirc templib_padcirc3 templib_padcirc2 templib_padcirc1 
                                      version mkdir netcdf netcdff XdmfCore XdmfUtils Xdmf
                                      ${MPI_Fortran_LIBRARIES})
    ELSE(NETCDF_WORKING AND NOT XDMF_WORKING)
        TARGET_LINK_LIBRARIES(padcirc templib_padcirc3 templib_padcirc2 templib_padcirc1 
                                      version mkdir ${MPI_Fortran_LIBRARIES})
    ENDIF(NETCDF_WORKING AND NOT XDMF_WORKING)
   
    ADD_DEPENDENCIES(padcirc          templib_padcirc3)
    ADD_DEPENDENCIES(templib_padcirc2 templib_padcirc1)
    ADD_DEPENDENCIES(templib_padcirc3 templib_padcirc2)
    ADD_DEPENDENCIES(templib_padcirc1 version)
    ADD_DEPENDENCIES(templib_padcirc1 mkdir)
    
    INSTALL(TARGETS padcirc RUNTIME DESTINATION bin)

ENDIF(BUILD_PADCIRC)
