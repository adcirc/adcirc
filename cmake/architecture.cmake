
#...Architecture Specifications. Determine if 
#   the system is 32 bit or 64 bit. If the 64 bit
#   integer pointer is not detected, an error is 
#   thrown during compile
SET(archdetect_c_code 
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
}"
)

FILE(WRITE "${CMAKE_BINARY_DIR}/CMakeFiles/architecture_check.c" "${archdetect_c_code}")

TRY_COMPILE(ARCH_TEST_COMPILE "${CMAKE_CURRENT_BINARY_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/architecture_check.c")
IF(ARCH_TEST_COMPILE)
    SET(ARCH 64)
ELSE(ARCH_TEST_COMPILE)
    SET(ARCH 32)
ENDIF(ARCH_TEST_COMPILE)
###########################################################################

INCLUDE(${CMAKE_SOURCE_DIR}/cmake/mpi_check.cmake)

###########################################################################
#...Compiler specific options
GET_FILENAME_COMPONENT(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
IF(Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  SET(Fortran_LINELENGTH_FLAG "-ffixed-line-length-none" CACHE STRING "Compiler specific flag to enable extended Fortran line length")

  # 64 bit array sizing
  IF(ARCH EQUAL 64)
    SET(Fortran_COMPILER_SPECIFIC_FLAG "-mcmodel=medium" CACHE STRING "Compiler specific flags")  
  ENDIF(ARCH EQUAL 64)

ELSEIF(Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort
  SET(Fortran_LINELENGTH_FLAG "-132" CACHE STRING "Compiler specific flag to enable extended Fortran line length")

  # Heap array allocation
  EXECUTE_PROCESS(COMMAND sh -c "ulimit -s" OUTPUT_VARIABLE STACKSIZE)
  STRING(STRIP "${STACKSIZE}" STACKSIZE_TRIMMED)
  MESSAGE(STATUS "The compiler flag -heap-arrays ${STACKSIZE_TRIMMED} is being added. This should also be used to compile netcdf-fortran.")
  
  # 64 bit array sizing
  IF(ARCH EQUAL 64)
    SET(Fortran_COMPILER_SPECIFIC_FLAG "-assume byterecl -mcmodel=medium -heap-arrays ${STACKSIZE_TRIMMED}" CACHE STRING "Compiler specific flags")  
  ELSE(ARCH EQUAL 64)
    SET(Fortran_COMPILER_SPECIFIC_FLAG "-assume byterecl -heap-arrays ${STACKSIZE_TRIMMED}" CACHE STRING "Compiler specific flags")
  ENDIF(ARCH EQUAL 64)

ELSEIF(Fortran_COMPILER_NAME MATCHES "pgf90.*")
  # pgf90
  SET(Fortran_LINELENGTH_FLAG "-Mextend" CACHE STRING "Compiler specific flag to enable extended Fortran line length")

  # 64 bit array sizing
  IF(ARCH EQUAL 64)
    SET(Fortran_COMPILER_SPECIFIC_FLAG "-Mlarge_arrays" CACHE STRING "Compiler specific flags")  
  ENDIF(ARCH EQUAL 64)

ELSE(Fortran_COMPILER_NAME MATCHES "gfortran.*")
  MESSAGE("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  MESSAGE("Fortran compiler: " ${Fortran_COMPILER_NAME})
  MESSAGE("No known predefined Fortran extended line length flag known. Please manually set the Fortran_LINELENGTH_FLAG")
  SET(Fortran_LINELENGTH_FLAG "" CACHE STRING "Compiler specific flag to enable extended Fortran line length")
  SET(Fortran_COMPILER_SPECIFIC_FLAG "" CACHE STRING "Compiler specific flags")
ENDIF(Fortran_COMPILER_NAME MATCHES "gfortran.*")
