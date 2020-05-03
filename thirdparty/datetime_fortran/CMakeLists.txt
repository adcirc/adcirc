project(datetime-fortran)
SET(CMAKE_POSITION_INDEPENDENT_CODE ON)
add_library(datetime ${CMAKE_CURRENT_SOURCE_DIR}/src/lib/datetime.f90 
                     ${CMAKE_CURRENT_SOURCE_DIR}/src/lib/mod_datetime.f90 
                     ${CMAKE_CURRENT_SOURCE_DIR}/src/lib/mod_timedelta.f90
                     ${CMAKE_CURRENT_SOURCE_DIR}/src/lib/mod_clock.f90 
                     ${CMAKE_CURRENT_SOURCE_DIR}/src/lib/mod_strftime.f90 
                     ${CMAKE_CURRENT_SOURCE_DIR}/src/lib/mod_constants.f90)
target_include_directories(datetime INTERFACE ${CMAKE_BINARY_DIR}/CMakeFiles/mod/datetime_fortran)
set_target_properties(datetime PROPERTIES
Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/CMakeFiles/mod/datetime_fortran
ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

if(WIN32)
  enable_language(CXX)
  add_library(strptime OBJECT strptime.cpp)
  set_target_properties(strptime PROPERTIES CXX_STANDARD 11 CXX_STANDARD_REQUIRED true)
  target_link_libraries(datetime PRIVATE strptime)
endif()

