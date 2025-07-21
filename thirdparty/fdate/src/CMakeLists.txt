include(GNUInstallDirs)
include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/SystemLink.cmake)

set(FDATE_SOURCES datetime_wrapper.cpp datetime.F90)

# Public headers that will be installed
set(FDATE_PUBLIC_HEADERS DateTime.hpp TimeDelta.hpp)

# Check for features
include(CheckCXXSourceCompiles)

# std::string_view
check_cxx_source_compiles("#include <string_view>
   int main() { std::string_view s; }" HAS_STRING_VIEW)

# Create interface library for common properties
add_library(fdate_interface INTERFACE)

# Create object library
add_library(fdate_objectlib OBJECT ${FDATE_SOURCES})
add_library(fdate $<TARGET_OBJECTS:fdate_objectlib>)

# Add alias for consistent naming
add_library(fdate::fdate ALIAS fdate)

# Link libraries
target_link_libraries(fdate_objectlib PUBLIC fdate::fdate_warnings
                                             fdate::fdate_options)

# Linker language
set_target_properties(fdate PROPERTIES LINKER_LANGUAGE Fortran)

# Compiler definitions
target_compile_definitions(
  fdate_objectlib PRIVATE $<$<BOOL:${HAS_STRING_VIEW}>:HAS_STRING_VIEW>
                          ONLY_C_LOCALE)

# Fortran module directory
set_target_properties(fdate_objectlib PROPERTIES Fortran_MODULE_DIRECTORY
                                                 ${CMAKE_BINARY_DIR}/mod)

# Set Fortran module directory for the main target too
set_target_properties(
  fdate
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod
             PUBLIC_HEADER "${FDATE_PUBLIC_HEADERS}"
             VERSION ${PROJECT_VERSION}
             SOVERSION ${PROJECT_VERSION_MAJOR})

# Include directories
target_include_directories(
  fdate_objectlib SYSTEM
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../thirdparty/date_hh)

target_include_directories(
  fdate PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
               $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# Installation
install(
  TARGETS fdate
  EXPORT fdateTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
          COMPONENT Runtime
          NAMELINK_COMPONENT Development
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT Development
  INCLUDES
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Install Fortran module files
install(
  DIRECTORY ${CMAKE_BINARY_DIR}/mod/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  COMPONENT Development
  FILES_MATCHING
  PATTERN "*.mod")

# Export targets
install(
  EXPORT fdateTargets
  FILE fdateTargets.cmake
  NAMESPACE fdate::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fdate
  COMPONENT Development)

# Create package config files
include(CMakePackageConfigHelpers)

# Generate the config file that includes the exports
configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/fdateConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/fdateConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fdate
  NO_SET_AND_CHECK_MACRO NO_CHECK_REQUIRED_COMPONENTS_MACRO)

# Generate the version file for the config file
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/fdateConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion)

# Install the configuration files
install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/fdateConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/fdateConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fdate
  COMPONENT Development)
