message( STATUS "System Architecture Detected: ${CMAKE_SYSTEM_PROCESSOR}")

include(${CMAKE_SOURCE_DIR}/cmake/mpi_check.cmake)

# ##############################################################################
# ...Compiler specific options
get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  # gfortran
  set(Fortran_LINELENGTH_FLAG
      "-ffixed-line-length-none"
      CACHE STRING
            "Compiler specific flag to enable extended Fortran line length")

  # 64 bit array sizing
  if(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64" )
    set(Fortran_COMPILER_SPECIFIC_FLAG
        "-mcmodel=medium"
        CACHE STRING "Compiler specific flags")
  else()
    set(Fortran_COMPILER_SPECIFIC_FLAG
        "" CACHE STRING "Compiler specific flag")
  endif()

elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel" OR ${CMAKE_Fortran_COMPILER_ID} STREQUAL "IntelLLVM")
  # ifort
  set(Fortran_LINELENGTH_FLAG
      "-132"
      CACHE STRING
            "Compiler specific flag to enable extended Fortran line length")

  # Heap array allocation
  execute_process(COMMAND sh -c "ulimit -s" OUTPUT_VARIABLE STACKSIZE)
  string(STRIP "${STACKSIZE}" STACKSIZE_TRIMMED)
  if("${STACKSIZE_TRIMMED}" STREQUAL "unlimited")
    message(STATUS "Unlimited stack size detected. No heap array flag added.")
    set(heaparray_FLAG " ")
  else("${STACKSIZE_TRIMMED}" STREQUAL "unlimited")
    message(
      STATUS
        "The compiler flag -heap-arrays ${STACKSIZE_TRIMMED} is being added. This should also be used to compile netcdf-fortran."
    )
    set(heaparray_FLAG "-heap-arrays ${STACKSIZE_TRIMMED}")
  endif("${STACKSIZE_TRIMMED}" STREQUAL "unlimited")

  # 64 bit array sizing
  if(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64")
    set(ifort_FLAG "${heaparray_FLAG} -assume byterecl -mcmodel=medium")
  else()
    set(ifort_FLAG "${heaparray_FLAG} -assume byterecl")
  endif()

  string(STRIP ${ifort_FLAG} ifort_FLAG_TRIMMED)
  set(Fortran_COMPILER_SPECIFIC_FLAG
      ${ifort_FLAG_TRIMMED}
      CACHE STRING "Compiler specific flags")

elseif(Fortran_COMPILER_NAME MATCHES "pgf90.*")
  # pgf90
  set(Fortran_LINELENGTH_FLAG
      "-Mextend"
      CACHE STRING
            "Compiler specific flag to enable extended Fortran line length")

  # 64 bit array sizing
  if(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64")
    set(Fortran_COMPILER_SPECIFIC_FLAG
        "-Mlarge_arrays"
        CACHE STRING "Compiler specific flags")
  endif()

else()
  message("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message(
    "No known predefined Fortran extended line length flag known. Please manually set the Fortran_LINELENGTH_FLAG"
  )
  set(Fortran_LINELENGTH_FLAG
      ""
      CACHE STRING
            "Compiler specific flag to enable extended Fortran line length")
  set(Fortran_COMPILER_SPECIFIC_FLAG
      ""
      CACHE STRING "Compiler specific flags")
endif()
