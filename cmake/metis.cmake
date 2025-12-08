# ######################################################################################################################
#
# ADCIRC - The ADvanced CIRCulation model Copyright (C) 1994-2025 R.A. Luettich, Jr., J.J. Westerink
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
add_library(adcirc_metis OBJECT ${METIS_SOURCES})
add_library(adcirc::metis ALIAS adcirc_metis)
target_include_directories(adcirc_metis PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/metis/Lib)
set_target_properties(adcirc_metis PROPERTIES EXCLUDE_FROM_ALL TRUE)

if(${CMAKE_C_COMPILER_ID} MATCHES "GNU")
  target_compile_options(
    adcirc_metis
    PRIVATE -Wno-implicit-function-declaration
            -Wno-incompatible-pointer-types
            -Wno-shift-op-parentheses
            -Wno-format-security)
elseif(${CMAKE_C_COMPILER_ID} MATCHES "AppleClang")
  target_compile_options(
    adcirc_metis
    PRIVATE -Wno-implicit-function-declaration
            -Wno-incompatible-pointer-types
            -Wno-shift-op-parentheses
            -Wno-format-security)
elseif(${CMAKE_C_COMPILER_ID} MATCHES "IntelLLVM")
  # For IntelLLVM, suppress warnings
  target_compile_options(adcirc_metis PRIVATE -Wno-incompatible-pointer-types -Wno-format-security
                                              -Wno-shift-op-parentheses)
elseif(${CMAKE_C_COMPILER_ID} MATCHES "Intel")
  # For Intel (classic), suppress warnings
  target_compile_options(adcirc_metis PRIVATE -diag-disable 167)
elseif(${CMAKE_C_COMPILER_ID} STREQUAL "PGI" OR ${CMAKE_C_COMPILER_ID} STREQUAL "NVHPC")
  # For PGI/NVHPC, suppress diagnostic warnings that are not part of ADCIRC
  target_compile_options(adcirc_metis PRIVATE --diag_suppress=167,177,550)
endif()
