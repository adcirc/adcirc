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

# Initialize lists for compile definitions and machine flags
set(ADCIRC_OPTION_FLAGS "")

if(ENABLE_POWELL)
  list(APPEND ADCIRC_OPTION_FLAGS "POWELL")
endif()

if(DEBUG_FULL_STACK)
  list(APPEND ADCIRC_OPTION_FLAGS "FULL_STACK")
endif()

if(DEBUG_ALL_TRACE)
  list(APPEND ADCIRC_OPTION_FLAGS "ALL_TRACE")
endif()

if(DEBUG_FLUSH_MESSAGES)
  list(APPEND ADCIRC_OPTION_FLAGS "FLUSH_MESSAGES")
endif()

if(DEBUG_LOG_LEVEL)
  list(APPEND ADCIRC_OPTION_FLAGS "EBUG")
endif()

if(DEBUG_GLOBALIO_TRACE)
  list(APPEND ADCIRC_OPTION_FLAGS "GLOBALIO_TRACE")
endif()

if(DEBUG_WRITER_TRACE)
  list(APPEND ADCIRC_OPTION_FLAGS "WRITER_TRACE")
endif()

if(DEBUG_WRITE_OUTPUT_TRACE)
  list(APPEND ADCIRC_OPTION_FLAGS "WRITE_OUTPUT_TRACE")
endif()

if(DEBUG_WIND_TRACE)
  list(APPEND ADCIRC_OPTION_FLAGS "WIND_TRACE")
endif()

if(DEBUG_WEIR_TRACE)
  list(APPEND ADCIRC_OPTION_FLAGS "WEIR_TRACE")
endif()

if(DEBUG_TVW_TRACE)
  list(APPEND ADCIRC_OPTION_FLAGS "TVW_TRACE")
endif()

if(DEBUG_VSMY_TRACE)
  list(APPEND ADCIRC_OPTION_FLAGS "VSMY_TRACE")
endif()

if(DEBUG_TIMESTEP_TRACE)
  list(APPEND ADCIRC_OPTION_FLAGS "TIMESTEP_TRACE")
endif()

if(DEBUG_SUBPREP_TRACE)
  list(APPEND ADCIRC_OPTION_FLAGS "SUBPREP_TRACE")
endif()

if(DEBUG_SUBDOMAIN_TRACE)
  list(APPEND ADCIRC_OPTION_FLAGS "SUBDOMAIN_TRACE")
endif()

if(DEBUG_READ_INPUT_TRACE)
  list(APPEND ADCIRC_OPTION_FLAGS "SUBDOMAIN_TRACE")
endif()

if(DEBUG_OWIWIND_TRACE)
  list(APPEND ADCIRC_OPTION_FLAGS "OWIWIND_TRACE")
endif()

if(DEBUG_NODALATTR_TRACE)
  list(APPEND ADCIRC_OPTION_FLAGS "NODALATTR_TRACE")
endif()

if(DEBUG_NETCDF_TRACE)
  list(APPEND ADCIRC_OPTION_FLAGS "NETCDF_TRACE")
endif()

if(DEBUG_MESSENGER_TRACE)
  list(APPEND ADCIRC_OPTION_FLAGS "MESSENGER_TRACE")
endif()

if(DEBUG_MESH_TRACE)
  list(APPEND ADCIRC_OPTION_FLAGS "MESH_TRACE")
endif()

if(DEBUG_HOTSTART_TRACE)
  list(APPEND ADCIRC_OPTION_FLAGS "HOTSTART_TRACE")
endif()

if(DEBUG_GLOBAL_TRACE)
  list(APPEND ADCIRC_OPTION_FLAGS "GLOBAL_TRACE")
endif()

if(DEBUG_HARM_TRACE)
  list(APPEND ADCIRC_OPTION_FLAGS "HARM_TRACE")
endif()

if(DEBUG_COLDSTART_TRACE)
  list(APPEND ADCIRC_OPTION_FLAGS "COLDSTART_TRACE")
endif()

if(DEBUG_COUPLE2SWAN_TRACE)
  list(APPEND ADCIRC_OPTION_FLAGS "COUPLE2SWAN_TRACE")
endif()

if(DEBUG_ADCIRC_TRACE)
  list(APPEND ADCIRC_OPTION_FLAGS "ADCIRC_TRACE")
endif()

if(DEBUG_HOLLAND)
  list(APPEND ADCIRC_OPTION_FLAGS "DEBUG_HOLLAND")
endif()

if(DEBUG_NWS14)
  list(APPEND ADCIRC_OPTION_FLAGS "DEBUG_NWS14")
endif()

if(DATETIME)
  list(APPEND ADCIRC_OPTION_FLAGS "DATETIME")
endif()

if(GRIB2)
  list(APPEND ADCIRC_OPTION_FLAGS "GRIB2API")
endif()

if(ADCIRC_NOF2008)
  list(APPEND ADCIRC_OPTION_FLAGS "NOF2008")
endif()

set(SWAN_FLAG "CSWAN")
