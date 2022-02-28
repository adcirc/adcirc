cmake_minimum_required( VERSION 3.9 )

list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )

project( gsw VERSION 3.05 LANGUAGES Fortran )

# Use solution folders in IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

### Use standard GNU installation directories
if ( NOT WIN32 )
    include( GNUInstallDirs )
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/gsw_version.f90.in ${CMAKE_CURRENT_BINARY_DIR}/gsw_version.f90)

### Configuration options
option( GSW_BUILD_STATIC_LIBS "Build static library" ON )
option( GSW_BUILD_SHARED_LIBS "Build shared library" OFF )

# Ensure at least one of GSW_BUILD_SHARED_LIBS and GSW_BUILD_STATIC_LIBS is set
if( NOT ( GSW_BUILD_STATIC_LIBS OR GSW_BUILD_SHARED_LIBS ) )
    message( STATUS "Niether GSW_BUILD_STATIC_LIBS nor GSW_BUILD_SHARED_LIBS is set.  Defaulting to GSW_BUILD_STATIC_LIBS=ON" )
    set( GSW_BUILD_STATIC_LIBS ON CACHE BOOL "Build static library" FORCE )
    set( GSW_BUILD_STATIC_LIBS ON )
endif()

message( STATUS "Option: GSW_BUILD_STATIC_LIBS: ${GSW_BUILD_STATIC_LIBS}" )
message( STATUS "Option: GSW_BUILD_SHARED_LIBS: ${GSW_BUILD_SHARED_LIBS}" )

# Set default build type to Release if not specified
if( NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE )
    message( STATUS "Setting default build type to 'Release'.  Set CMAKE_BUILD_TYPE variable to change build types." )
    set_property( CACHE CMAKE_BUILD_TYPE PROPERTY VALUE "Release" )
endif()

# Postfix for debug libraries so they can live alongside release libraries
set( CMAKE_DEBUG_POSTFIX ".debug" CACHE STRING "Debug file extension" )

## Package compiler flags
include( compiler_flags )

### Global compilation properties
set( INCLUDE_DIR ${PROJECT_NAME} ) # path relative to <prefix>/include/ to install headers/modules
set( CMAKE_Fortran_MODULE_DIRECTORY ${INCLUDE_DIR} )
set( CMAKE_POSITION_INDEPENDENT_CODE ON )

### Define Library Targets
add_subdirectory(modules)
add_subdirectory(toolbox)
list( APPEND SRCS ${modules_srcs} ${toolbox_srcs} )

# Use a common object library for building shared and static targets
add_library( ${PROJECT_NAME}_objects OBJECT ${SRCS} )
target_compile_definitions( ${PROJECT_NAME}_objects PUBLIC ${PUBLIC_FLAGS}
                                                    PRIVATE ${PRIVATE_FLAGS} )
set_property(TARGET ${PROJECT_NAME}_objects PROPERTY FOLDER gsw)

# Add static lib target
if( GSW_BUILD_STATIC_LIBS )
    add_library( ${PROJECT_NAME}_static STATIC ${CMAKE_CURRENT_BINARY_DIR}/gsw_version.f90 $<TARGET_OBJECTS:${PROJECT_NAME}_objects> )
    list( APPEND LIB_TARGETS ${PROJECT_NAME}_static )
endif()

# Add shared lib target
if( GSW_BUILD_SHARED_LIBS )
    add_library( ${PROJECT_NAME}_shared SHARED ${CMAKE_CURRENT_BINARY_DIR}/gsw_version.f90 $<TARGET_OBJECTS:${PROJECT_NAME}_objects> )
    list( APPEND LIB_TARGETS ${PROJECT_NAME}_shared )
endif()

# Set common lib target properties
set_target_properties( ${LIB_TARGETS} PROPERTIES OUTPUT_NAME ${PROJECT_NAME} )
foreach( _tgt IN LISTS LIB_TARGETS )
    target_compile_definitions( ${_tgt} PUBLIC "${PUBLIC_FLAGS}" )
    target_include_directories( ${_tgt}
        PUBLIC
            $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/${INCLUDE_DIR}> )
    set_property(TARGET ${_tgt} PROPERTY FOLDER gsw)
endforeach()

# Enable ctest framework and add tests
enable_testing()
if ( NOT WIN32 )
    add_subdirectory( test EXCLUDE_FROM_ALL )
else()
    add_subdirectory( test )
endif()

### Install
install( TARGETS ${LIB_TARGETS} EXPORT ${PROJECT_NAME}Exports
         ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
         INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${INCLUDE_DIR} )
install( DIRECTORY ${PROJECT_BINARY_DIR}/${INCLUDE_DIR}/
         DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${INCLUDE_DIR} )

### Package config
include( CMakePackageConfigHelpers )

export( EXPORT ${PROJECT_NAME}Exports NAMESPACE ${PROJECT_NAME}:: FILE ${PROJECT_NAME}-targets.cmake )

set( CONFIG_INSTALL_DESTINATION share/${PROJECT_NAME}/cmake )

# gsw-config.cmake
configure_package_config_file( cmake/PackageConfig.cmake.in ${PROJECT_NAME}-config.cmake
                               INSTALL_DESTINATION ${CONFIG_INSTALL_DESTINATION} )
install( FILES ${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config.cmake
         DESTINATION ${CONFIG_INSTALL_DESTINATION} )

# gsw-config-version.cmake
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion )
install( FILES ${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake
         DESTINATION ${CONFIG_INSTALL_DESTINATION} )

# gsw-targets.cmake and gsw-targets-<build-type>.cmake
install( EXPORT ${PROJECT_NAME}Exports NAMESPACE ${PROJECT_NAME}::
         FILE ${PROJECT_NAME}-targets.cmake
         DESTINATION ${CONFIG_INSTALL_DESTINATION} )
