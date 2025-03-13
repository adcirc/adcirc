# Adds strict warning flags to *.F90 files in the ADCIRC repository. These flags help diagnose potential issues
macro(add_strict_compiler_flags)
  if(${ADCIRC_STRICT_COMPILER_FLAGS})

    if(${CMAKE_Fortran_COMPILER_ID} MATCHES "IntelLLVM")
      set(STRICT_FLAGS "-warn all -diag-enable remark -implicit-none")
    elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
      set(STRICT_FLAGS
          "-Werror -Wall -Wextra -Wconversion -pedantic -fimplicit-none -Wuninitialized -Wsurprising -Wuse-without-only -Wimplicit-procedure -Winteger-division -Wconversion-extra"
      )
    else()
      message(STATUS "No strict compiler flags defined for ${CMAKE_Fortran_COMPILER_ID}. No action taken.")
    endif()

    foreach(SOURCE_FILE ${ARGN})
      if(${SOURCE_FILE} MATCHES ".*\\.F90$")
        set_source_files_properties(${SOURCE_FILE} PROPERTIES COMPILE_FLAGS "${STRICT_FLAGS}")
      endif()
    endforeach()

  endif()

endmacro()
