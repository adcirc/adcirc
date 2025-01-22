add_library(kdtree ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/KDTREE2/kdtree2.F)
set(LOCAL_COMPILER_FLAGS "${Fortran_LINELENGTH_FLAG} ${Fortran_COMPILER_SPECIFIC_FLAG}")
string(STRIP ${LOCAL_COMPILER_FLAGS} LOCAL_COMPILER_FLAGS)
set_target_properties(kdtree PROPERTIES COMPILE_FLAGS ${LOCAL_COMPILER_FLAGS})
set_target_properties(kdtree PROPERTIES EXCLUDE_FROM_ALL TRUE)
set_target_properties(kdtree PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/mod/kdtree)
set_target_properties(kdtree PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
