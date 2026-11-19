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
# ...Lenient date parsing used by all ADCIRC executables (depends only on fdate)
add_library(adcirc_date_parsing OBJECT ${CMAKE_CURRENT_SOURCE_DIR}/src/date_parsing.F90)
add_library(adcirc::date_parsing ALIAS adcirc_date_parsing)
adcirc_set_module_directory(adcirc_date_parsing)
target_include_directories(adcirc_date_parsing PRIVATE ${CMAKE_BINARY_DIR}/mod)
add_dependencies(adcirc_date_parsing fdate)
set_target_properties(adcirc_date_parsing PROPERTIES EXCLUDE_FROM_ALL TRUE)
