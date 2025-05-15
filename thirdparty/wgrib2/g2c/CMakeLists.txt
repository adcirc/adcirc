# This is the main CMake file for NCEPLIBS-g2c.
#
# Mark Potts, Kyle Gerheiser, Ed Hartnett, Eric Engle
cmake_minimum_required(VERSION 3.15)

# Read the VERSION file.
file(STRINGS "VERSION" pVersion)

# Set up the project.
project(g2c VERSION ${pVersion} LANGUAGES C)

# Provide install directories according to GNU standards.
include(GNUInstallDirs)

# Find CMake code we need.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

# Handle user options.
option(ENABLE_DOCS "Enable generation of doxygen-based documentation." OFF)
option(USE_PNG "Use PNG library"      ON)
option(USE_Jasper "Use Jasper library"   ON)
option(USE_OpenJPEG "Use OpenJPEG library" OFF)
option(USE_AEC "Use LibAEC library" OFF)
option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(BUILD_STATIC_LIBS "Build static libraries" ON)
option(FTP_TEST_FILES "Fetch and test with files on FTP site." OFF)
option(FTP_LARGE_TEST_FILES "Fetch and test with very large files on FTP site." OFF)
option(FTP_EXTRA_TEST_FILES "Test with more large files fetched via FTP." OFF)
option(LOGGING "Turn on internal logging messages. Only useful to g2c developers." OFF)
option(PTHREADS "Turn on thread-safty with pthreads." OFF)
option(UTILS "Build and install some utility programs." ON)
option(BUILD_G2C "Build the g2c file-based API." ON)

# Developers can use this option to specify a local directory which
# holds the test files. They will be copied instead of fetching the
# files via FTP.
SET(TEST_FILE_DIR "." CACHE STRING "Check this directory for test files before using FTP.")
message(STATUS "Finding test data files in directory ${TEST_FILE_DIR}.")

# Set this to better handle files > 2 GB.
add_compile_definitions(_LARGEFILE64_SOURCE)

# Set pre-processor symbol if logging is desired.
if(LOGGING)
  add_definitions(-DLOGGING)
endif()

# Set pre-processor symbol if thread safety is desired.
if(PTHREADS)
  add_definitions(-DPTHREADS)
endif()

# The user must select either Jasper of OpenJPEG
if(USE_Jasper AND USE_OpenJPEG)
  message(FATAL_ERROR "Either Jasper or OpenJPEG should be used, not both.")
endif()

# Check the build type.
if(NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|Release|RelWithDebInfo|MinSizeRel)$")
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
                                               "MinSizeRel" "RelWithDebInfo")
endif()

# Find required packages to use PNG.
if(USE_PNG)
  find_package(ZLIB REQUIRED)
  find_package(PkgConfig REQUIRED)
  find_package(PNG REQUIRED)
  set(PNG_WIDTH_MAX 100000000 CACHE STRING "Default width of PNG grid.")
  set(PNG_HEIGHT_MAX 100000 CACHE STRING "Default height of PNG grid.")
  message(STATUS "\tWill build with PNG grid width = ${PNG_WIDTH_MAX} and height = ${PNG_HEIGHT_MAX}")
else()
  message(STATUS "Will not build PNG support")
endif()

# Find Jasper if it is desired.
if(USE_Jasper)
  find_package(Jasper 2.0.25 REQUIRED)
  set(JASPER_MAX_MEM 268435456 CACHE STRING "Default size of Jasper memory buffer.")
  message(STATUS "\tWill build with Jasper memory buffer size = ${JASPER_MAX_MEM}")
else()
  message(STATUS "Will not build Jasper support")
endif()

# Find OpenJpeg if it is desired.
if(USE_OpenJPEG)
  find_package(OpenJPEG REQUIRED)
else()
  message(STATUS "Will not build OpenJPEG support")
endif()

# Turn on this pre-processor symbol to get JPEG code and testing.
if(USE_Jasper OR USE_OpenJPEG)
  add_definitions(-DJPEG)
endif()

# Find required packages to use AEC.
if(USE_AEC)
  find_package(libaec 1.0.6 REQUIRED)
  add_definitions(-DAEC)
else()
  message(STATUS "Will not build AEC support")
endif()

# Set the compiler flags.
if(CMAKE_C_COMPILER_ID MATCHES "^(Intel|IntelLLVM)$")
  set(CMAKE_C_FLAGS "-g -traceback ${CMAKE_C_FLAGS}")
  set(CMAKE_C_FLAGS_RELEASE "-O3")
  set(CMAKE_C_FLAGS_DEBUG "-O0")
elseif(CMAKE_C_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
  set(CMAKE_C_FLAGS "-g ${CMAKE_C_FLAGS}")
  set(CMAKE_C_FLAGS_RELEASE "-O3")
  set(CMAKE_C_FLAGS_DEBUG "-O0 -ggdb -Wall")
endif()

# Set this to better handle files > 2 GB.
add_compile_definitions(_LARGEFILE64_SOURCE)

set(lib_name ${PROJECT_NAME})

# Build the code in the source directory.
add_subdirectory(src)

# Build the code in the utils directory.
if(BUILD_G2C)
  if(UTILS)
    add_subdirectory(utils)
  endif()
endif()

# Determine whether or not to generate documentation.
if(ENABLE_DOCS)
  find_package(Doxygen REQUIRED)
endif()
add_subdirectory(docs)

# Run unit tests.
include(CTest)
if(BUILD_TESTING)
  add_subdirectory(tests)
endif()
