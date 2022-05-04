# ...Architecture Specifications. Determine if the system is 32 bit or 64 bit.
# If the 64 bit integer pointer is not detected, an error is thrown during
# compile
set(archdetect_c_code
        "#include <stdint.h>
int main()
{
#if INTPTR_MAX == INT64_MAX
        return 0;
#elif INTPTR_MAX == INT32_MAX
#error 32-bit max integer pointer
#else
#error Unknown pointer size or missing size macros
#endif
}")

file(WRITE "${CMAKE_BINARY_DIR}/CMakeFiles/architecture_check.c"
        "${archdetect_c_code}")

try_compile(ARCH_TEST_COMPILE "${CMAKE_CURRENT_BINARY_DIR}"
        "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/architecture_check.c")
if (ARCH_TEST_COMPILE)
  set(ARCH 64)
else (ARCH_TEST_COMPILE)
  set(ARCH 32)
endif (ARCH_TEST_COMPILE)
# ##############################################################################

include(${CMAKE_SOURCE_DIR}/cmake/mpi_check.cmake)

# ##############################################################################
# ...Compiler specific options
get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  set(Fortran_LINELENGTH_FLAG
          "-ffixed-line-length-none"
          CACHE STRING
          "Compiler specific flag to enable extended Fortran line length")

  # 64 bit array sizing
  if (ARCH EQUAL 64)
    set(Fortran_COMPILER_SPECIFIC_FLAG
            "-mcmodel=medium"
            CACHE STRING "Compiler specific flags")
  endif (ARCH EQUAL 64)

elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort
  set(Fortran_LINELENGTH_FLAG
          "-132"
          CACHE STRING
          "Compiler specific flag to enable extended Fortran line length")

  # Heap array allocation
  execute_process(COMMAND sh -c "ulimit -s" OUTPUT_VARIABLE STACKSIZE)
  string(STRIP "${STACKSIZE}" STACKSIZE_TRIMMED)
  if ("${STACKSIZE_TRIMMED}" STREQUAL "unlimited")
    message(STATUS "Unlimited stack size detected. No heap array flag added.")
    set(heaparray_FLAG " ")
  else ("${STACKSIZE_TRIMMED}" STREQUAL "unlimited")
    message(
            STATUS
            "The compiler flag -heap-arrays ${STACKSIZE_TRIMMED} is being added. This should also be used to compile netcdf-fortran."
    )
    set(heaparray_FLAG "-heap-arrays ${STACKSIZE_TRIMMED}")
  endif ("${STACKSIZE_TRIMMED}" STREQUAL "unlimited")

  # 64 bit array sizing
  if (ARCH EQUAL 64)
    set(ifort_FLAG "${heaparray_FLAG} -assume byterecl -mcmodel=medium")
  else (ARCH EQUAL 64)
    set(ifort_FLAG "${heaparray_FLAG} -assume byterecl")
  endif (ARCH EQUAL 64)

  string(STRIP ${ifort_FLAG} ifort_FLAG_TRIMMED)
  set(Fortran_COMPILER_SPECIFIC_FLAG
          ${ifort_FLAG_TRIMMED}
          CACHE STRING "Compiler specific flags")

elseif (Fortran_COMPILER_NAME MATCHES "pgf90.*")
  # pgf90
  set(Fortran_LINELENGTH_FLAG
          "-Mextend"
          CACHE STRING
          "Compiler specific flag to enable extended Fortran line length")

  # 64 bit array sizing
  if (ARCH EQUAL 64)
    set(Fortran_COMPILER_SPECIFIC_FLAG
            "-Mlarge_arrays"
            CACHE STRING "Compiler specific flags")
  endif (ARCH EQUAL 64)

else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
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
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")
