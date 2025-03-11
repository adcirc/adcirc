macro(adcirc_enable_ipo)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT result OUTPUT output)
  if(result)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
  else()
    message(WARNING "IPO is not supported: ${output}")
  endif()
endmacro()

option(ADCIRC_ENABLE_IPO "Enable Interprocedural Optimization (IPO)" OFF)
if(ADCIRC_ENABLE_IPO)
  adcirc_enable_ipo()
endif()
mark_as_advanced(ADCIRC_ENABLE_IPO)
