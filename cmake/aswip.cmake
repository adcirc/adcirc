
IF(BUILD_ASWIP)


    SET(ASWIP_SOURCES ${CMAKE_SOURCE_DIR}/src/sizes.F 
                      ${CMAKE_SOURCE_DIR}/src/global.F 
                      ${CMAKE_SOURCE_DIR}/src/global_3dvs.F 
                      ${CMAKE_SOURCE_DIR}/src/boundaries.F 
                      ${CMAKE_SOURCE_DIR}/src/hashtable.F 
                      ${CMAKE_SOURCE_DIR}/src/nodalattr.F 
                      ${CMAKE_SOURCE_DIR}/src/mesh.F 
                      ${CMAKE_SOURCE_DIR}/src/wind.F 
                      ${CMAKE_SOURCE_DIR}/src/owiwind.F 
                      ${CMAKE_SOURCE_DIR}/thirdparty/KDTREE2/kdtree2.F 
                      ${CMAKE_SOURCE_DIR}/src/owi_ice.F
                      ${CMAKE_SOURCE_DIR}/wind/vortex.F 
                      ${CMAKE_SOURCE_DIR}/wind/aswip.F )

    ADD_EXECUTABLE(aswip ${ASWIP_SOURCES})

    addCompilerFlags(aswip ${ADDITIONAL_FLAGS_ASWIP})
    addLibVersion(aswip)

    INSTALL(TARGETS aswip RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

ENDIF(BUILD_ASWIP)
