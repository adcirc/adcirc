
IF(BUILD_ASWIP)
    
    SET(ASWIP_SOURCES src/sizes.F src/global.F src/global_3dvs.F src/boundaries.F src/nodalattr.F 
                      src/mesh.F src/wind.F src/owiwind.F KDTREE2/kdtree2.F src/owi_ice.F 
                      wind/vortex.F wind/aswip_1.0.3.F )

    ADD_EXECUTABLE(aswip ${ASWIP_SOURCES})

    addCompilerFlags(aswip ${ADDITIONAL_FLAGS_ASWIP})
    addLibVersion(aswip) 

    INSTALL(TARGETS aswip RUNTIME DESTINATION bin)

ENDIF(BUILD_ASWIP)
