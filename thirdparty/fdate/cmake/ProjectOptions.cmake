include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/SystemLink.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/LibFuzzer.cmake)
include(CMakeDependentOption)
include(CheckCXXCompilerFlag)

macro(fdate_supports_sanitizers)

  # If we are on Darwin arm64, we can't use sanitizers
  if(CMAKE_SYSTEM_NAME STREQUAL "Darwin" AND CMAKE_SYSTEM_PROCESSOR STREQUAL
                                             "arm64")
    set(SUPPORTS_ASAN OFF)
    set(SUPPORTS_UBSAN OFF)
  else()
    if((CMAKE_CXX_COMPILER_ID MATCHES ".*Clang.*" OR CMAKE_CXX_COMPILER_ID
                                                     MATCHES ".*GNU.*")
       AND NOT WIN32)
      set(SUPPORTS_UBSAN ON)
    else()
      set(SUPPORTS_UBSAN OFF)
    endif()

    if((CMAKE_CXX_COMPILER_ID MATCHES ".*Clang.*" OR CMAKE_CXX_COMPILER_ID
                                                     MATCHES ".*GNU.*")
       AND WIN32)
      set(SUPPORTS_ASAN OFF)
    else()
      set(SUPPORTS_ASAN ON)
    endif()
  endif()
endmacro()

macro(fdate_setup_options)

  fdate_supports_sanitizers()

  option(FDATE_ENABLE_COVERAGE "Enable coverage reporting" OFF)

  if(NOT PROJECT_IS_TOP_LEVEL OR NOT FDATE_DEVELOPER_MODE)
    option(fdate_DEVELOPER_MODE "Enable developer mode" OFF)
    option(FDATE_ENABLE_HARDENING "Enable hardening" OFF)
    cmake_dependent_option(
      FDATE_ENABLE_GLOBAL_HARDENING
      "Attempt to push hardening options to built dependencies" OFF
      FDATE_ENABLE_HARDENING OFF)
    option(FDATE_ENABLE_IPO "Enable IPO/LTO" OFF)
    option(FDATE_WARNINGS_AS_ERRORS "Treat Warnings As Errors" OFF)
    option(FDATE_ENABLE_USER_LINKER "Enable user-selected linker" OFF)
    option(FDATE_ENABLE_SANITIZER_ADDRESS "Enable address sanitizer" OFF)
    option(FDATE_ENABLE_SANITIZER_LEAK "Enable leak sanitizer" OFF)
    option(FDATE_ENABLE_SANITIZER_UNDEFINED "Enable undefined sanitizer" OFF)
    option(FDATE_ENABLE_SANITIZER_THREAD "Enable thread sanitizer" OFF)
    option(FDATE_ENABLE_SANITIZER_MEMORY "Enable memory sanitizer" OFF)
    option(FDATE_ENABLE_UNITY_BUILD "Enable unity builds" OFF)
    # option(FDATE_ENABLE_CLANG_TIDY "Enable clang-tidy" OFF)
    option(FDATE_ENABLE_CPPCHECK "Enable cpp-check analysis" OFF)
    option(FDATE_ENABLE_PCH "Enable precompiled headers" OFF)
    option(FDATE_ENABLE_CACHE "Enable ccache" OFF)
  else()
    option(FDATE_DEVELOPER_MODE "Enable developer mode" ON)
    option(FDATE_ENABLE_HARDENING "Enable hardening" ON)
    cmake_dependent_option(
      FDATE_ENABLE_GLOBAL_HARDENING
      "Attempt to push hardening options to built dependencies" ON
      FDATE_ENABLE_HARDENING OFF)
    option(FDATE_ENABLE_IPO "Enable IPO/LTO" OFF)
    option(FDATE_WARNINGS_AS_ERRORS "Treat Warnings As Errors" ON)
    option(FDATE_ENABLE_USER_LINKER "Enable user-selected linker" OFF)
    option(FDATE_ENABLE_SANITIZER_ADDRESS "Enable address sanitizer"
           ${SUPPORTS_ASAN})
    option(FDATE_ENABLE_SANITIZER_LEAK "Enable leak sanitizer" OFF)
    option(FDATE_ENABLE_SANITIZER_UNDEFINED "Enable undefined sanitizer"
           ${SUPPORTS_UBSAN})
    option(FDATE_ENABLE_SANITIZER_THREAD "Enable thread sanitizer" OFF)
    option(FDATE_ENABLE_SANITIZER_MEMORY "Enable memory sanitizer" OFF)
    option(FDATE_ENABLE_UNITY_BUILD "Enable unity builds" OFF)
    # option(FDATE_ENABLE_CLANG_TIDY "Enable clang-tidy" ON)
    option(FDATE_ENABLE_CPPCHECK "Enable cpp-check analysis" ON)
    option(FDATE_ENABLE_PCH "Enable precompiled headers" OFF)
    option(FDATE_ENABLE_CACHE "Enable ccache" ON)
  endif()

  if(NOT PROJECT_IS_TOP_LEVEL OR NOT fdate_DEVELOPER_MODE)
    mark_as_advanced(
      FDATE_DEVELOPER_MODE
      FDATE_ENABLE_HARDENING
      FDATE_ENABLE_GLOBAL_HARDENING
      FDATE_ENABLE_IPO
      FDATE_WARNINGS_AS_ERRORS
      FDATE_ENABLE_USER_LINKER
      FDATE_ENABLE_SANITIZER_ADDRESS
      FDATE_ENABLE_SANITIZER_LEAK
      FDATE_ENABLE_SANITIZER_UNDEFINED
      FDATE_ENABLE_SANITIZER_THREAD
      FDATE_ENABLE_SANITIZER_MEMORY
      FDATE_ENABLE_UNITY_BUILD
      # FDATE_ENABLE_CLANG_TIDY
      FDATE_ENABLE_CPPCHECK
      FDATE_ENABLE_COVERAGE
      FDATE_ENABLE_PCH
      FDATE_ENABLE_CACHE)
  endif()

  if(FDATE_DEVELOPER_MODE)
    fdate_check_libfuzzer_support(LIBFUZZER_SUPPORTED)
    if(LIBFUZZER_SUPPORTED
       AND (FDATE_ENABLE_SANITIZER_ADDRESS
            OR FDATE_ENABLE_SANITIZER_THREAD
            OR FDATE_ENABLE_SANITIZER_UNDEFINED))
      set(DEFAULT_FUZZER ON)
    else()
      set(DEFAULT_FUZZER OFF)
    endif()

    option(FDATE_BUILD_FUZZ_TESTS "Enable fuzz testing executable"
           ${DEFAULT_FUZZER})
  endif()

endmacro()

macro(fdate_global_options)
  if(FDATE_ENABLE_IPO)
    include(cmake/InterproceduralOptimization.cmake)
    fdate_enable_ipo()
  endif()

  fdate_supports_sanitizers()

  if(FDATE_ENABLE_HARDENING AND FDATE_ENABLE_GLOBAL_HARDENING)
    include(cmake/Hardening.cmake)
    if(NOT SUPPORTS_UBSAN
       OR FDATE_ENABLE_SANITIZER_UNDEFINED
       OR FDATE_ENABLE_SANITIZER_ADDRESS
       OR FDATE_ENABLE_SANITIZER_THREAD
       OR FDATE_ENABLE_SANITIZER_LEAK)
      set(ENABLE_UBSAN_MINIMAL_RUNTIME FALSE)
    else()
      set(ENABLE_UBSAN_MINIMAL_RUNTIME TRUE)
    endif()
    # message("${FDATE_ENABLE_HARDENING} ${ENABLE_UBSAN_MINIMAL_RUNTIME}
    # ${FDATE_ENABLE_SANITIZER_UNDEFINED}")
    fdate_enable_hardening(fdate_options ON ${ENABLE_UBSAN_MINIMAL_RUNTIME})
  endif()
endmacro()

macro(fdate_local_options)
  if(PROJECT_IS_TOP_LEVEL)
    include(cmake/StandardProjectSettings.cmake)
  endif()

  add_library(fdate_warnings INTERFACE)
  add_library(fdate_options INTERFACE)

  include(cmake/CompilerWarnings.cmake)
  fdate_set_project_warnings(fdate_warnings ${FDATE_WARNINGS_AS_ERRORS} "" ""
                             "" "")

  if(FDATE_ENABLE_USER_LINKER)
    include(cmake/Linker.cmake)
    fdate_configure_linker(fdate_options)
  endif()

  include(cmake/Sanitizers.cmake)
  fdate_enable_sanitizers(
    fdate_options ${FDATE_ENABLE_SANITIZER_ADDRESS}
    ${FDATE_ENABLE_SANITIZER_LEAK} ${FDATE_ENABLE_SANITIZER_UNDEFINED}
    ${FDATE_ENABLE_SANITIZER_THREAD} ${FDATE_ENABLE_SANITIZER_MEMORY})

  set_target_properties(fdate_options PROPERTIES UNITY_BUILD
                                                 ${FDATE_ENABLE_UNITY_BUILD})

  if(FDATE_ENABLE_PCH)
    target_precompile_headers(fdate_options INTERFACE <vector> <string>
                              <utility>)
  endif()

  if(FDATE_ENABLE_CACHE)
    include(cmake/Cache.cmake)
    fdate_enable_cache()
  endif()

  include(cmake/StaticAnalyzers.cmake)
  # if(FDATE_ENABLE_CLANG_TIDY) fdate_enable_clang_tidy(fdate_options
  # ${fdate_WARNINGS_AS_ERRORS}) endif()

  if(FDATE_ENABLE_CPPCHECK)
    fdate_enable_cppcheck(
      ${FDATE_WARNINGS_AS_ERRORS} "" # override cppcheck options
    )
  endif()

  if(FDATE_ENABLE_COVERAGE)
    include(cmake/Tests.cmake)
    fdate_enable_coverage(fdate_options)
  endif()

  if(fdate_WARNINGS_AS_ERRORS)
    check_cxx_compiler_flag("-Wl,--fatal-warnings" LINKER_FATAL_WARNINGS)
    if(LINKER_FATAL_WARNINGS)
      # This is not working consistently, so disabling for now
      # target_link_options(fdate_options INTERFACE -Wl,--fatal-warnings)
    endif()
  endif()

  if(FDATE_ENABLE_HARDENING AND NOT FDATE_ENABLE_GLOBAL_HARDENING)
    include(cmake/Hardening.cmake)
    if(NOT SUPPORTS_UBSAN
       OR FDATE_ENABLE_SANITIZER_UNDEFINED
       OR FDATE_ENABLE_SANITIZER_ADDRESS
       OR FDATE_ENABLE_SANITIZER_THREAD
       OR FDATE_ENABLE_SANITIZER_LEAK)
      set(ENABLE_UBSAN_MINIMAL_RUNTIME FALSE)
    else()
      set(ENABLE_UBSAN_MINIMAL_RUNTIME TRUE)
    endif()
    fdate_enable_hardening(fdate_options OFF ${ENABLE_UBSAN_MINIMAL_RUNTIME})
  endif()

endmacro()
