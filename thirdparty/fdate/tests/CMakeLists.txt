# ##############################################################################
# ...CXX TESTS
# ##############################################################################
file(GLOB_RECURSE TestList ${CMAKE_CURRENT_SOURCE_DIR}/cxx/TEST_*.cpp)

foreach(test ${TestList})
  get_filename_component(test_name ${test} NAME_WE)
  add_executable(${test_name} ${test})
  add_test(
    NAME ${test_name}
    COMMAND ${test_name}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  target_include_directories(${test_name}
                             PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../src)

  target_link_libraries(${test_name} PRIVATE fdate::fdate_warnings
                                             fdate::fdate_options)

  if(FDATE_ENABLE_COVERAGE)
    set_target_properties(
      ${test_name} PROPERTIES COMPILE_FLAGS ${FDATE_COVERAGE_COMPILE_FLAGS}
                              LINK_FLAGS ${FDATE_COVERAGE_LINK_FLAGS})
  endif()

  target_include_directories(
    ${test_name} SYSTEM
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../thirdparty/date_hh)

  target_link_libraries(${test_name} PRIVATE Catch2::Catch2WithMain)
  set_target_properties(
    ${test_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                            ${CMAKE_CURRENT_BINARY_DIR}/tests)
  add_dependencies(${test_name} Catch2::Catch2WithMain)
endforeach()
# ##############################################################################

# ##############################################################################
# Check CXX results against Python results
# ##############################################################################
find_package(Python3 COMPONENTS Interpreter)
if(Python3_FOUND)
  add_executable(datetime_interface
                 ${CMAKE_CURRENT_SOURCE_DIR}/cxx/cxx_test_interface.cpp)

  add_executable(datetime_ftn_interface
                 ${CMAKE_CURRENT_SOURCE_DIR}/fortran/ftn_test_interface.F90)

  target_include_directories(datetime_interface
                             PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../src)
  target_include_directories(
    datetime_interface SYSTEM
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../thirdparty/date_hh)

  target_include_directories(datetime_ftn_interface
                             PRIVATE ${CMAKE_BINARY_DIR}/mod)

  target_link_libraries(datetime_interface PRIVATE fdate fdate::fdate_warnings
                                                   fdate::fdate_options)

  target_link_libraries(
    datetime_ftn_interface PRIVATE fdate fdate::fdate_warnings
                                   fdate::fdate_options)

  set_target_properties(datetime_ftn_interface PROPERTIES LINKER_LANGUAGE
                                                          Fortran)

  set(FDATE_FUZZ_TEST_DURATION 10)

  # Add a test which calls the Python script that runs random tests against
  # python's datetime library
  add_test(
    NAME TEST_cxx_python_fuzz_parse
    COMMAND
      ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/python/python_fuzz_wrapper.py --timeout
      ${FDATE_FUZZ_TEST_DURATION} --test parse --exe
      ${CMAKE_BINARY_DIR}/tests/datetime_interface)

  add_test(
    NAME TEST_cxx_python_fuzz_arithmetic
    COMMAND
      ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/python/python_fuzz_wrapper.py --timeout
      ${FDATE_FUZZ_TEST_DURATION} --test arithmetic --exe
      ${CMAKE_BINARY_DIR}/tests/datetime_interface)

  add_test(
    NAME TEST_ftn_python_fuzz_parse
    COMMAND
      ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/python/python_fuzz_wrapper.py --timeout
      ${FDATE_FUZZ_TEST_DURATION} --test parse --exe
      ${CMAKE_BINARY_DIR}/tests/datetime_ftn_interface)

  add_test(
    NAME TEST_ftn_python_fuzz_arithmetic
    COMMAND
      ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/python/python_fuzz_wrapper.py --timeout
      ${FDATE_FUZZ_TEST_DURATION} --test arithmetic --exe
      ${CMAKE_BINARY_DIR}/tests/datetime_ftn_interface)

endif()

# ##############################################################################
# ... Fortran TESTS
# ##############################################################################
file(GLOB_RECURSE TestList ${CMAKE_CURRENT_SOURCE_DIR}/fortran/TEST_ftn_*.F90)
foreach(test ${TestList})
  get_filename_component(test_name ${test} NAME_WE)
  add_executable(${test_name} ${test})
  add_test(
    NAME ${test_name}
    COMMAND ${test_name}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  target_include_directories(${test_name} PRIVATE ${CMAKE_BINARY_DIR}/mod)

  target_link_libraries(${test_name} PRIVATE fdate fdate::fdate_warnings
                                             fdate::fdate_options)

  # Set the linker language to Fortran
  set_target_properties(${test_name} PROPERTIES LINKER_LANGUAGE Fortran)

  set_target_properties(
    ${test_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                            ${CMAKE_CURRENT_BINARY_DIR}/tests)
  add_dependencies(${test_name} fdate)
endforeach()
# ##############################################################################
