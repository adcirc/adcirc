AC_INIT([adcirc],[v56.0.2],[zcobell@thewaterinstitute.org])
AC_CONFIG_HEADERS(config.h)
AC_CONFIG_SRCDIR([src])
AM_INIT_AUTOMAKE([subdir-objects])
AC_CONFIG_MACRO_DIR([m4])
LT_INIT()

dnl Check for base compilers
AC_PROG_CC
AC_PROG_F77 dnl Note that we'll end up always using the f90 compiler
AC_PROG_FC
AC_FC_LINE_LENGTH([132])

dnl Add option for --with-swan
ac_with_swan=no
AC_ARG_WITH(swan, AS_HELP_STRING([--with-swan],[Build coupled ADCIRC+SWAN]),[
if test "x$with_swan" = "xyes" ; then
    ac_with_swan=yes
else
    ac_with_swan=no
fi
])

dnl Add option for --with-mpi
ac_have_mpi=no
AC_ARG_WITH(mpi, AS_HELP_STRING([--with-mpi],[Build with parallel support]),[
if test "x$with_mpi" = "xyes" ; then
    ac_build_mpi=yes
else
    ac_build_mpi=no
fi
])
  

dnl Add option for --with-netcdf
AC_ARG_WITH(netcdf,
  AS_HELP_STRING([--with-netcdf],
                 [Use netCDF. Default: auto-detect]), [
case "$with_netcdf" in
  yes|no)
    : # Nothing special to do here
    ;;
  *)
    if test ! -d "$withval" ; then
      AC_MSG_ERROR([--with-netcdf path does not point to a directory])
    fi
    LIBNETCDF_PATH="-L$with_netcdf/lib"
    CFLAGS="-I$with_netcdf/include $CFLAGS"
    FCFLAGS="-I$with_netcdf/include $FCFLAGS"
  esac
])

dnl If we are using netcdf, check the configuration here
nc_has_fortran=0
if test "x$with_netcdf" != "xno" ; then

  if test "x$with_netcdf" == "x" || test "x$with_netcdf" == "xyes" ; then
    LIBNETCDF_PATH="/usr"
    CFLAGS="-I/usr/include $CFLAGS"
    FCFLAGS="-I/usr/include $FCFLAGS"
  fi

  AC_LANG_PUSH([C])
  AC_CHECK_HEADERS(netcdf.h)
  AC_CHECK_LIB([netcdf],[nc_open],[LIBNETCDF_LIBS="$LIBNETCDF_PATH -lnetcdf"])
  AC_LANG_POP()

  AC_LANG_PUSH([Fortran])
  AC_MSG_CHECKING([for netcdf-fortran])
  AC_COMPILE_IFELSE([
    AC_LANG_PROGRAM(,[[
          use netcdf            
          implicit none
                    ]])],
   [nc_has_fortran=yes],[nc_has_fortran=no]
  )
  AC_MSG_RESULT([$nc_has_fortran])

  if test "x$nc_has_fortran" != "xyes" ; then
      AC_MSG_ERROR([netcdf-fortran is required to enable netcdf. Try configuring with --with-netcdf=no])
  fi

  AC_LANG_POP()

fi

dnl If we are using MPI, then check for the compiler wrapper
if test "x$ac_build_mpi" = "xyes" ; then
    AC_LANG_PUSH([Fortran])
    AX_MPI([ac_have_mpi=yes],[ac_have_mpi=no])
    AC_LANG_POP()
fi

dnl If building swan, we need perl
AS_IF([test x$ac_with_swan = "xyes"],
    [
     AC_CHECK_PROG([HAS_PERL],[perl],[yes],[no])
     AS_IF([test x$HAS_PERL != xyes],AC_MSG_ERROR([Perl required to build SWAN source code]),[])
    ],[])

dnl To generate the version file, we need python
AC_CHECK_PROG([HAS_PYTHON],[python],[yes],[no])
AS_IF([test x$HAS_PYTHON != xyes],AC_MSG_ERROR([Python require to generate the version information]),[])

dnl Set the conditional values in all Makefile.am files 
AM_CONDITIONAL([BUILD_MPI], [test x$ac_have_mpi = xyes])
AM_CONDITIONAL([BUILD_NETCDF], [test x$nc_has_fortran = xyes])
AM_CONDITIONAL([BUILD_SWAN], [test x$ac_with_swan])

dnl Generate the SWAN source files via the preprocessor
AS_IF([test x$ac_with_swan = "xyes"],
    [
     AC_CONFIG_COMMANDS([adcswan], 
        [
         perl thirdparty/swan/switch.pl -unix -adcirc thirdparty/swan/*.ftn thirdparty/swan/*.ftn90
         mkdir -p autotools/lib/libswan_serial_source
         mv thirdparty/swan/*.f thirdparty/swan/*.f90 autotools/lib/libswan_serial_source/.
        ]
     )
    ],[]
)
AS_IF([test x$ac_with_swan = "xyes" && test x$ac_have_mpi = xyes],
    [
     AC_CONFIG_COMMANDS([padcswan], 
        [
         perl thirdparty/swan/switch.pl -unix -adcirc -pun thirdparty/swan/*.ftn thirdparty/swan/*.ftn90
         mkdir -p autotools/lib/libswan_parallel_source
         mv thirdparty/swan/*.f thirdparty/swan/*.f90 autotools/lib/libswan_parallel_source/.
        ]
     )
    ],[]
)

AC_CONFIG_COMMANDS([version],
  [scripts/adcirc_version.py --create-version-file --directory $ac_top_srcdir]
)

dnl Set the Makefiles to be created
AC_CONFIG_FILES([
  Makefile
  autotools/Makefile
  autotools/lib/Makefile
  autotools/lib/libswan_serial1/Makefile
  autotools/lib/libadcswan/Makefile
  autotools/lib/libswan_serial2/Makefile
  autotools/lib/metis/Makefile
  autotools/bin/Makefile
  autotools/bin/adcirc/Makefile
  autotools/bin/adcprep/Makefile
  autotools/bin/padcirc/Makefile
  autotools/bin/adcswan/Makefile
])

dnl Write the output
AC_OUTPUT
