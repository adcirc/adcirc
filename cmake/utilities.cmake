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
if(BUILD_UTILITIES)

  add_executable(adccmp util/adccmp.F)
  add_executable(p15 wind/p15.F)
  add_executable(owi22 wind/owi22.F)
  add_executable(build13 util/build13.F)
  add_executable(buildstwave23 util/buildstwave23.F)
  add_executable(hot2asc util/hot2asc.F)
  add_executable(inflate util/inflate.F)
  add_executable(hstime util/hstime.F)

  add_executable(adcircResultsComparison util/adcircResultsComparison.F90)
  addcompilerflags(adcircResultsComparison)

  addcompilerflags(adccmp)
  addcompilerflags(p15)
  addcompilerflags(owi22)
  addcompilerflags(build13)
  addcompilerflags(buildstwave23)
  addcompilerflags(hot2asc)
  addcompilerflags(inflate)
  addcompilerflags(hstime)

  addnetcdflibraries(hstime)
  addnetcdflibraries(adcircResultsComparison)

  # Some of these utilities are very old and no longer updated. We will
  # pass some compiler flags to suppress warnings. If this happened in the
  # main code, we would fix the code instead of suppressing the warning.
  if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    # Disable warnings for deleted features
    target_compile_options(p15 PRIVATE -std=legacy)
  endif()


  install(
    TARGETS adccmp
            p15
            owi22
            build13
            buildstwave23
            hot2asc
            inflate
            hstime
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

endif(BUILD_UTILITIES)
