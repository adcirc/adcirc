add_library(mkdir ${CMAKE_CURRENT_SOURCE_DIR}/prep/mkdir.c)
target_include_directories(mkdir PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
set_target_properties(mkdir PROPERTIES EXCLUDE_FROM_ALL TRUE)
