# ######################################################################################################################
#
# ADCIRC - The ADvanced CIRCulation model Copyright (C) 1994-2023 R.A. Luettich, Jr., J.J. Westerink
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General
# Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along with this program.  If not, see
# <http://www.gnu.org/licenses/>.
#
# ######################################################################################################################
set(METIS_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/coarsen.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/fm.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/initpart.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/match.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/ccgraph.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/memory.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/pmetis.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/pqueue.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/refine.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/util.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/timing.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/debug.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/bucketsort.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/graph.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/stat.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/kmetis.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/kwayrefine.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/kwayfm.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/balance.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/ometis.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/srefine.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/sfm.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/separator.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/mincover.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/mmd.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/mesh.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/meshpart.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/frename.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/fortran.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/myqsort.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/compress.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/parmetis.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/estmem.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/mpmetis.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/mcoarsen.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/mmatch.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/minitpart.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/mbalance.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/mrefine.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/mutil.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/mfm.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/mkmetis.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/mkwayrefine.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/mkwayfmh.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/mrefine2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/minitpart2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/mbalance2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/mfm2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/kvmetis.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/kwayvolrefine.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/kwayvolfm.c
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/metis/Lib/subdomains.c)
add_library(metis STATIC ${METIS_SOURCES})
target_include_directories(metis PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/metis/Lib)
set_target_properties(metis PROPERTIES EXCLUDE_FROM_ALL TRUE)

if(${CMAKE_C_COMPILER_ID} MATCHES "GNU")
  # For GCC >= 10, we need to add -Wno-implicit-function-declaration and -Wno-incompatible-pointer-types
  if(${CMAKE_C_COMPILER_VERSION} VERSION_GREATER 10 OR ${CMAKE_C_COMPILER_VERSION} VERSION_EQUAL 10)
    set(ADDITIONAL_METIS_COMPILER_FLAGS
        "${ADDITIONAL_METIS_COMPILER_FLAGS} -Wno-implicit-function-declaration -Wno-incompatible-pointer-types -Wno-shift-op-parentheses -Wno-format-security"
    )
    message(STATUS "Adding additional compiler flags to metis: ${ADDITIONAL_METIS_COMPILER_FLAGS}")
    set_target_properties(metis PROPERTIES COMPILE_FLAGS ${ADDITIONAL_METIS_COMPILER_FLAGS})
  endif()
elseif(${CMAKE_C_COMPILER_ID} MATCHES "IntelLLVM")
  # For IntelLLVM, we need to add -Wno-incompatible-pointer-types -Wno-format-security -Wno-shift-op-parentheses
  set(ADDITIONAL_METIS_COMPILER_FLAGS
      "${ADDITIONAL_METIS_COMPILER_FLAGS} -Wno-incompatible-pointer-types -Wno-format-security -Wno-shift-op-parentheses"
  )
  message(STATUS "Adding additional compiler flags to metis: ${ADDITIONAL_METIS_COMPILER_FLAGS}")
  set_target_properties(metis PROPERTIES COMPILE_FLAGS ${ADDITIONAL_METIS_COMPILER_FLAGS})
elseif(${CMAKE_C_COMPILER_ID} MATCHES "Intel")
  # For Intel (classic), we need to add -Wno-incompatible-pointer-types -Wno-format-security -Wno-shift-op-parentheses
  set(ADDITIONAL_METIS_COMPILER_FLAGS "${ADDITIONAL_METIS_COMPILER_FLAGS} -diag-disable 167")
  message(STATUS "Adding additional compiler flags to metis: ${ADDITIONAL_METIS_COMPILER_FLAGS}")
  set_target_properties(metis PROPERTIES COMPILE_FLAGS ${ADDITIONAL_METIS_COMPILER_FLAGS})
elseif(${CMAKE_C_COMPILER_ID} MATCHES "NVHPC")
  # For nvc, suppress diagnostic warnings that are not part of ADCIRC
  set(ADDITIONAL_METIS_COMPILER_FLAGS
      "${ADDITIONAL_METIS_COMPILER_FLAGS} --display_error_number --diag_suppress 550 --diag_suppress 177 --diag_suppress 167"
  )
  set_target_properties(metis PROPERTIES COMPILE_FLAGS ${ADDITIONAL_METIS_COMPILER_FLAGS})
endif()
