# This is the main CMake file for the src directory of NCEPLIBS-g2c.
#
# Ed Hartnett 12/29/22
# Eric Engle 1/10/23

set(lib_name ${PROJECT_NAME})

# Specify the source code files for the library.
set(src_files grib2_int.h cmplxpack.c compack.c
comunpack.c drstemplates.c g2_addfield.c g2_addgrid.c
g2_addlocal.c g2_create.c g2_free.c g2_getfld.c
g2_gribend.c g2_info.c g2_miss.c g2_unpack1.c
g2_unpack2.c g2_unpack3.c g2_unpack4.c g2_unpack5.c
g2_unpack6.c g2_unpack7.c gbits.c getdim.c
getpoly.c gridtemplates.c int_power.c misspack.c
mkieee.c pack_gp.c pdstemplates.c rdieee.c
reduce.c seekgb.c simpack.c simunpack.c specpack.c
specunpack.c util.c)

if(BUILD_G2C)
  set(src_files ${src_files}  g2cutil.c g2ccsv.c g2cparams.c g2cfile.c
    g2cdegrib2.c g2cprod.c g2cinq.c g2cindex.c g2cio.c g2ccompare.c)
endif()

# Create object library.  The object library will be used to create
# static and shared libraries.
add_library(${lib_name}_objlib OBJECT ${src_files})
set_property(TARGET ${lib_name}_objlib PROPERTY C_STANDARD 11)
target_include_directories(${lib_name}_objlib PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
if(BUILD_SHARED_LIBS)
  set_property(TARGET ${lib_name}_objlib PROPERTY POSITION_INDEPENDENT_CODE 1)
endif()
add_library(${PROJECT_NAME}::${lib_name} ALIAS ${lib_name}_objlib)

# Create static object library from object library.
if(BUILD_STATIC_LIBS)
  add_library(${lib_name}_static STATIC $<TARGET_OBJECTS:${lib_name}_objlib>)
  set_target_properties(${lib_name}_static PROPERTIES OUTPUT_NAME ${lib_name})
  set_target_properties(${lib_name}_static PROPERTIES EXPORT_NAME ${lib_name})
  set_target_properties(${lib_name}_static PROPERTIES EXPORT_NAME ${lib_name}_static)
  list(APPEND LIB_TARGETS ${lib_name}_static)
endif()

# Create shared object library from object library.
if(BUILD_SHARED_LIBS)
  add_library(${lib_name}_shared SHARED $<TARGET_OBJECTS:${lib_name}_objlib>)
  set_target_properties(${lib_name}_shared PROPERTIES OUTPUT_NAME ${lib_name})
  set_target_properties(${lib_name}_shared PROPERTIES EXPORT_NAME ${lib_name})
  set_target_properties(${lib_name}_shared PROPERTIES SOVERSION 0)
  list(APPEND LIB_TARGETS ${lib_name}_shared)
endif()

# Build with PNG.
if(PNG_FOUND)
  message(STATUS "PNG include directories ${PNG_INCLUDE_DIR}")
  message(STATUS "PNG libraries ${PNG_LIBRARIES}")
  target_sources(${lib_name}_objlib PRIVATE decenc_png.c
                                     pngpack.c
                                     pngunpack.c)
  target_compile_definitions(${lib_name}_objlib PRIVATE USE_PNG)
  target_link_libraries(${lib_name}_objlib PRIVATE ${PNG_LIBRARIES})
  target_include_directories(${lib_name}_objlib PRIVATE "${PNG_INCLUDE_DIRS}")
  set(G2_PNG_ENABLED ON)
  if(BUILD_SHARED_LIBS)
    target_link_libraries(${lib_name}_shared PRIVATE ${PNG_LIBRARIES})
    target_include_directories(${lib_name}_shared PRIVATE "${PNG_INCLUDE_DIRS}")
  endif()
  if(BUILD_STATIC_LIBS)
    target_link_libraries(${lib_name}_static PRIVATE ${PNG_LIBRARIES})
    target_include_directories(${lib_name}_static PRIVATE "${PNG_INCLUDE_DIRS}")
  endif()
endif()

# Build with Jasper.
if(JASPER_FOUND)
  message(STATUS "JASPER include directories ${JASPER_INCLUDE_DIR}")
  message(STATUS "JASPER libraries ${JASPER_LIBRARIES} version ${JASPER_VERSION}")
  target_sources(${lib_name}_objlib PRIVATE decenc_jpeg2000.c
                                     jpcpack.c
                                     jpcunpack.c)
  target_compile_definitions(${lib_name}_objlib PRIVATE USE_JPEG2000)
  target_link_libraries(${lib_name}_objlib PRIVATE ${JASPER_LIBRARIES})
  target_include_directories(${lib_name}_objlib PRIVATE "${JASPER_INCLUDE_DIR}")
  set(G2_JPEG2000_ENABLED ON)
  if(JASPER_VERSION_STRING VERSION_GREATER_EQUAL "3")
    message(STATUS "      JASPER is version 3.0 or greater")
    add_definitions(-DJASPER3)
  endif()
  if(BUILD_SHARED_LIBS)
    target_link_libraries(${lib_name}_shared PRIVATE ${JASPER_LIBRARIES})
    target_include_directories(${lib_name}_shared PRIVATE "${JASPER_INCLUDE_DIRS}")
  endif()
  if(BUILD_STATIC_LIBS)
    target_link_libraries(${lib_name}_static PRIVATE ${JASPER_LIBRARIES})
    target_include_directories(${lib_name}_static PRIVATE "${JASPER_INCLUDE_DIRS}")
  endif()
endif()

# Build with OpenJPEG.
if(OPENJPEG_FOUND)
  message(STATUS "OpenJPEG include directories ${OpenJPEG_INCLUDE_DIRS}")
  message(STATUS "OpenJPEG libraries ${OpenJPEG_LIBRARIES}")
  target_sources(${lib_name}_objlib PRIVATE decenc_openjpeg.c
                                     jpcpack.c
                                     jpcunpack.c)
  target_compile_definitions(${lib_name}_objlib PRIVATE USE_OPENJPEG)
  target_link_libraries(${lib_name}_objlib PRIVATE ${OpenJPEG_LIBRARIES})
  target_include_directories(${lib_name}_objlib PRIVATE "${OpenJPEG_INCLUDE_DIRS}")
  set(G2_JPEG2000_ENABLED ON)
  if(BUILD_SHARED_LIBS)
    target_link_libraries(${lib_name}_shared PRIVATE ${OpenJPEG_LIBRARIES})
    target_include_directories(${lib_name}_shared PRIVATE "${OpenJPEG_INCLUDE_DIRS}")
  endif()
  if(BUILD_STATIC_LIBS)
    target_link_libraries(${lib_name}_static PRIVATE ${OpenJPEG_LIBRARIES})
    target_include_directories(${lib_name}_static PRIVATE "${OpenJPEG_INCLUDE_DIRS}")
  endif()
endif()

# Build with AEC.
if(LIBAEC_FOUND)
  set(AEC_LIBRARIES ${LIBAEC_LIBRARIES})
  set(AEC_INCLUDE_DIR ${LIBAEC_INCLUDE_DIRS})
  message(STATUS "AEC include directories ${AEC_INCLUDE_DIR}")
  message(STATUS "AEC libraries ${AEC_LIBRARIES}")
  target_sources(${lib_name}_objlib PRIVATE decenc_aec.c
                                            aecpack.c
                                            aecunpack.c)
  target_compile_definitions(${lib_name}_objlib PRIVATE USE_AEC)
  target_link_libraries(${lib_name}_objlib PRIVATE ${AEC_LIBRARIES})
  target_include_directories(${lib_name}_objlib PRIVATE "${AEC_INCLUDE_DIR}")
  set(G2_AEC_ENABLED ON)
  if(BUILD_SHARED_LIBS)
    target_link_libraries(${lib_name}_shared PRIVATE ${AEC_LIBRARIES})
    target_include_directories(${lib_name}_shared PRIVATE "${AEC_INCLUDE_DIR}")
  endif()
  if(BUILD_STATIC_LIBS)
    target_link_libraries(${lib_name}_static PRIVATE ${AEC_LIBRARIES})
    target_include_directories(${lib_name}_static PRIVATE "${AEC_INCLUDE_DIR}")
  endif()
endif()

# Build with pthreads.
if(PTHREADS)
  target_link_libraries(${lib_name}_objlib PRIVATE pthread)
  if(BUILD_SHARED_LIBS)
    target_link_libraries(${lib_name}_shared PRIVATE pthread)
  endif()
  if(BUILD_STATIC_LIBS)
    target_link_libraries(${lib_name}_static PRIVATE pthread)
  endif()
endif()

# Create this header file, which has the contents of the VERSION file
# substituted in.
CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/grib2.h.in" "${CMAKE_CURRENT_BINARY_DIR}/grib2.h"
  @ONLY)

target_include_directories(${lib_name}_objlib PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR};${CMAKE_CURRENT_SOURCE_DIR}>"
                                                     $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>)
if(BUILD_SHARED_LIBS)
  target_include_directories(${lib_name}_shared PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR};${CMAKE_CURRENT_SOURCE_DIR}>"
                                                       $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>)
endif()
if(BUILD_STATIC_LIBS)
  target_include_directories(${lib_name}_static PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR};${CMAKE_CURRENT_SOURCE_DIR}>"
                                                       $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>)
endif()

if (NOT WIN32)
  target_link_libraries(${lib_name}_objlib INTERFACE m)
  if(BUILD_SHARED_LIBS)
    target_link_libraries(${lib_name}_shared INTERFACE m)
  endif()
  if(BUILD_STATIC_LIBS)
    target_link_libraries(${lib_name}_static INTERFACE m)
  endif()
endif()

if (WIN32)
  target_link_libraries(${lib_name}_objlib INTERFACE ws2_32)
endif()

# Install g2c library
install(
  TARGETS ${LIB_TARGETS}
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Install generated grib2.h header
install(FILES ${CMAKE_BINARY_DIR}/src/grib2.h DESTINATION ${CMAKE_INSTALL_PREFIX}/include)

# Package config
include(CMakePackageConfigHelpers)
set(CONFIG_INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})

export(EXPORT ${PROJECT_NAME}Exports
  NAMESPACE ${PROJECT_NAME}::
  FILE ${PROJECT_NAME}-targets.cmake)

# g2c-config.cmake
configure_package_config_file(
  ${PROJECT_SOURCE_DIR}/cmake/PackageConfig.cmake.in ${CMAKE_BINARY_DIR}/${PROJECT_NAME}-config.cmake
  INSTALL_DESTINATION ${CONFIG_INSTALL_DESTINATION})
install(FILES ${CMAKE_BINARY_DIR}/${PROJECT_NAME}-config.cmake
  DESTINATION ${CONFIG_INSTALL_DESTINATION})

# g2c-config-version.cmake
write_basic_package_version_file(
  ${CMAKE_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion)
install(FILES ${CMAKE_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake
  DESTINATION ${CONFIG_INSTALL_DESTINATION})

# g2c-targets.cmake and g2c-targets-<build-type>.cmake
install(EXPORT ${PROJECT_NAME}Exports
  NAMESPACE ${PROJECT_NAME}::
  FILE ${PROJECT_NAME}-targets.cmake
  DESTINATION ${CONFIG_INSTALL_DESTINATION})

# CSV data files.
if(BUILD_G2C)
  set(CSV_FILES Template.txt CodeFlag.txt)
  install(FILES ${CSV_FILES}
    DESTINATION ${CMAKE_INSTALL_PREFIX}/share/g2c/)
endif()
