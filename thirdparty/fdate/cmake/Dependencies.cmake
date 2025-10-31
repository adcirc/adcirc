function(mark_as_advanced_wildcard variable_prefix)
  get_cmake_property(CACHE_VARS VARIABLES)
  foreach(CACHE_VAR ${CACHE_VARS})
    if(CACHE_VAR MATCHES "^${variable_prefix}")
      mark_as_advanced(${CACHE_VAR})
    endif()
  endforeach()
endfunction()

macro(fdate_setup_dependencies)

  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/CPM.cmake)

  # ############################################################################
  # Catch2
  # ############################################################################
  if(FDATE_ENABLE_TESTING)

    # Mark the CPM variables as advanced
    mark_as_advanced_wildcard("CPM_")

    cpmaddpackage(
      NAME
      Catch2
      GITHUB_REPOSITORY
      catchorg/Catch2
      VERSION
      3.7.0
      EXCLUDE_FROM_ALL
      SYSTEM)
    mark_as_advanced_wildcard("CATCH_")

    # Disable cppcheck and clang-tidy on Catch2
    set_target_properties(Catch2 PROPERTIES CXX_CPPCHECK "" CXX_CLANG_TIDY "")
    set_target_properties(Catch2WithMain PROPERTIES CXX_CPPCHECK ""
                                                    CXX_CLANG_TIDY "")

  endif()

  # Mark the FetchContent variables as advanced
  mark_as_advanced_wildcard("FETCHCONTENT_")
  mark_as_advanced("pugixml_DIR")

endmacro()
