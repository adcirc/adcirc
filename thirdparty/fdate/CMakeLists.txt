cmake_minimum_required(VERSION 3.14)
project(fdate VERSION 0.0.1)

# Packaging and uninstall support only make sense when fdate is built on its own; when it is
# embedded with add_subdirectory they would leak into the parent project.
if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(FDATE_IS_TOP_LEVEL ON)
else()
  set(FDATE_IS_TOP_LEVEL OFF)
endif()
include(CheckCXXCompilerFlag)

check_cxx_compiler_flag("-std=c++17" COMPILER_SUPPORTS_CXX17)
if(COMPILER_SUPPORTS_CXX17)
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
else()
  message(
    FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} does not support C++17.")
endif()

enable_language(CXX)
enable_language(Fortran)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Dependencies.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/PreventInSourceBuilds.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/ProjectOptions.cmake)
fdate_setup_options()
fdate_global_options()
fdate_local_options()
add_library(fdate::fdate_options ALIAS fdate_options)
add_library(fdate::fdate_warnings ALIAS fdate_warnings)
fdate_setup_dependencies()

add_subdirectory(src)

option(FDATE_ENABLE_TESTING "Enable FDATE testing" OFF)
if(FDATE_ENABLE_TESTING)
  enable_testing()
  add_subdirectory(tests)
endif()

if(FDATE_IS_TOP_LEVEL)
  # Configure CPack for packaging
  set(CPACK_PACKAGE_NAME "fdate")
  set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY ${PROJECT_DESCRIPTION})
  set(CPACK_PACKAGE_VENDOR "Zach Cobell")
  set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
  set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
  include(CPack)

  # Add uninstall target
  if(NOT TARGET uninstall)
    configure_file(
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
      "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake" IMMEDIATE @ONLY)

    add_custom_target(
      uninstall COMMAND ${CMAKE_COMMAND} -P
                        ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
  endif()
endif()
