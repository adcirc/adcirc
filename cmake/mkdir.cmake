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
add_library(adcirc_mkdir OBJECT ${CMAKE_CURRENT_SOURCE_DIR}/prep/mkdir.c)
add_library(adcirc::mkdir ALIAS adcirc_mkdir)
target_include_directories(adcirc_mkdir PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
set_target_properties(adcirc_mkdir PROPERTIES EXCLUDE_FROM_ALL TRUE)
