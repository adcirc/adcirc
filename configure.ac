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

dnl Set the conditional values in all Makefile.am files 
AM_CONDITIONAL([HAVE_MPI], [test x$ac_have_mpi = xyes])
AM_CONDITIONAL([HAVE_NETCDF], [test x$nc_has_fortran = xyes])

dnl Set the Makefiles to be created
AC_CONFIG_FILES([
  Makefile
  autotools/Makefile
  autotools/adcirc/Makefile
  autotools/metis/Makefile
  autotools/adcprep/Makefile
  autotools/padcirc/Makefile
])

dnl Write the output
AC_OUTPUT
