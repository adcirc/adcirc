if(BUILD_UTILITIES)

  add_executable(adccmp util/adccmp.F)
  add_executable(p15 wind/p15.F)
  add_executable(owi22 wind/owi22.F)
  add_executable(build13 util/build13.F)
  add_executable(buildstwave23 util/buildstwave23.F)
  add_executable(hot2asc util/hot2asc.F)
  add_executable(inflate util/inflate.F)
  add_executable(hstime util/hstime.F)

  if(NETCDF_WORKING)
    add_executable(adcircResultsComparison util/adcircResultsComparison.F90)
    addcompilerflags(adcircResultsComparison)
    addnetcdf(adcircResultsComparison)
  endif(NETCDF_WORKING)

  addcompilerflags(adccmp)
  addcompilerflags(p15)
  addcompilerflags(owi22)
  addcompilerflags(build13)
  addcompilerflags(buildstwave23)
  addcompilerflags(hot2asc)
  addcompilerflags(inflate)
  addcompilerflags(hstime)

  install(
    TARGETS adccmp
            p15
            owi22
            build13
            buildstwave23
            hot2asc
            inflate
            hstime
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

endif(BUILD_UTILITIES)
