add_library(version ${CMAKE_BINARY_DIR}/CMakeFiles/version_cmake.F)
if(WIN32)
  if(EXISTS ${CMAKE_SOURCE_DIR}/version.F)
    add_custom_command(
      OUTPUT ${CMAKE_BINARY_DIR}/CMakeFiles/version_cmake.F
      COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/version.F
              ${CMAKE_BINARY_DIR}/CMakeFiles/version_cmake.F
      COMMENT "Generating ADCIRC version...")
  elseif(EXISTS ${CMAKE_SOURCE_DIR}/version_default.F)
    add_custom_command(
      OUTPUT ${CMAKE_BINARY_DIR}/CMakeFiles/version_cmake.F
      COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/version_default.F
              ${CMAKE_BINARY_DIR}/CMakeFiles/version_cmake.F
      COMMENT "Generating ADCIRC version...")
  endif(EXISTS ${CMAKE_SOURCE_DIR}/version.F)
else(WIN32)
  add_custom_command(
    OUTPUT ${CMAKE_BINARY_DIR}/CMakeFiles/version_cmake.F
    COMMAND ./adcircVersion.sh ${CMAKE_SOURCE_DIR} >/dev/null
    COMMAND cp ${CMAKE_SOURCE_DIR}/version.F
            ${CMAKE_BINARY_DIR}/CMakeFiles/version_cmake.F
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/scripts
    COMMENT "Generating ADCIRC version...")
endif(WIN32)
set_target_properties(version PROPERTIES Fortran_MODULE_DIRECTORY
                                         CMakeFiles/version_mod)
set_target_properties(
  version
  PROPERTIES COMPILE_FLAGS
             "${Fortran_LINELENGTH_FLAG} ${Fortran_COMPILER_SPECIFIC_FLAG}")
set_target_properties(version PROPERTIES EXCLUDE_FROM_ALL TRUE)
