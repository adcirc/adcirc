
ADD_LIBRARY(version version.F)
#ADD_LIBRARY(version ${CMAKE_BINARY_DIR}/version.F)
#ADD_CUSTOM_COMMAND( OUTPUT ${CMAKE_BINARY_DIR}/version_cmake.F
#    COMMAND ./generateVersion.sh ${CMAKE_BINARY_DIR}/version_cmake.F
#    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/cmake
#    COMMENT "Generating ADCIRC version...")
SET_TARGET_PROPERTIES(version PROPERTIES Fortran_MODULE_DIRECTORY CMakeFiles/version_mod)
