#
#  Find OpenJPEG
#
#  OpenJPEG_FOUND        - if false, do not try to use OpenJPEG
#  OpenJPEG_INCLUDE_DIRS - the OpenJPEG include directories
#  OpenJPEG_LIBRARIES    - the libraries to link against to use OpenJPEG
#
#  This file is part of NCEPLIBS-g2c. Distributed under the LGPL v3.0.

find_package(PkgConfig QUIET)
pkg_check_modules(OpenJPEG_PKGCONF QUIET libopenjp2)

set(OpenJPEG_VERSION ${OpenJPEG_PKGCONF_VERSION})

find_path(OpenJPEG_INCLUDE_DIRS
  NAMES openjpeg.h
  HINTS ${OpenJPEG_PKGCONF_INCLUDEDIR} ${OpenJPEG_PKGCONF_INCLUDE_DIRS}
)

find_library(OpenJPEG_LIBRARIES
  NAMES openjp2
  HINTS ${OpenJPEG_PKGCONF_LIBDIR} ${OpenJPEG_PKGCONF_LIBRARY_DIRS}
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OpenJPEG
    FOUND_VAR OpenJPEG_FOUND
    REQUIRED_VARS OpenJPEG_LIBRARIES OpenJPEG_INCLUDE_DIRS
    VERSION_VAR OpenJPEG_VERSION
)

mark_as_advanced(OpenJPEG_INCLUDE_DIRS OpenJPEG_LIBRARIES)
