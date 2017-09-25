
IF(BUILD_UTILITIES)
    
    ADD_EXECUTABLE(adccmp util/adccmp.F)
    ADD_EXECUTABLE(p15 wind/p15.F)
    ADD_EXECUTABLE(owi22 wind/owi22.F)
    ADD_EXECUTABLE(build13 util/build13.F)
    ADD_EXECUTABLE(buildstwave23 util/buildstwave23.F)
    ADD_EXECUTABLE(hot2asc util/hot2asc.F)
    ADD_EXECUTABLE(inflate util/inflate.F)
    ADD_EXECUTABLE(hstime util/hstime.F)

    IF(NETCDF_WORKING)
        ADD_EXECUTABLE(adcircResultCompare util/adcircResultCompare.F90)
        addCompilerFlags(adcircResultCompare)
        addNetCDF(adcircResultCompare)
    ENDIF(NETCDF_WORKING)

    addCompilerFlags(adccmp)
    addCompilerFlags(p15)
    addCompilerFlags(owi22)
    addCompilerFlags(build13)
    addCompilerFlags(buildstwave23)
    addCompilerFlags(hot2asc)
    addCompilerFlags(inflate)
    addCompilerFlags(hstime)
    
    INSTALL(TARGETS adccmp p15 owi22 build13 buildstwave23 hot2asc inflate hstime 
            RUNTIME DESTINATION bin)

ENDIF(BUILD_UTILITIES)
