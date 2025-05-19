AC_INIT([adcirc],[v56.0.2],[zcobell@thewaterinstitute.org])

# Find out about the host we're building on.
AC_CANONICAL_HOST

# Find out about the target we're building for.
AC_CANONICAL_TARGET

AM_INIT_AUTOMAKE([subdir-objects foreign dist-zip])
AM_MAINTAINER_MODE
AC_CONFIG_HEADERS(config.h)
AC_CONFIG_SRCDIR([src])
AM_SILENT_RULES([yes])
AC_CONFIG_MACRO_DIR([m4])
LT_INIT()

dnl Check for base compilers
AC_PROG_CC
AC_PROG_FC
AC_PROG_F77($FC) 
AC_FC_LINE_LENGTH([132])

dnl gfortran-10 errors out if argument mismatch (SWAN)
dnl See if Fotran compiler supports -fallow-argument-mismatch for legacy code
AC_LANG_PUSH([Fortran])
AX_CHECK_COMPILE_FLAG([-fallow-argument-mismatch],
                      [FFLAGS="${FFLAGS} -fallow-argument-mismatch" FCFLAGS="${FCFLAGS} -fallow-argument-mismatch"],
                      [])
AC_LANG_POP([Fortran])
AC_LANG_PUSH([C])
AX_CHECK_COMPILE_FLAG([-Wno-incompatble-pointer-types],
                      [CFLAGS="${CFLAGS} -Wno-incompatible-pointer-types"],
                      [])
AC_LANG_POP([C])

dnl Add option for --enable-swan
ac_enable_swan=no
AC_ARG_ENABLE(swan, AS_HELP_STRING([--enable-swan],[Build coupled ADCIRC+SWAN]),[
if test "x$enable_swan" = "xyes" ; then
    ac_enable_swan=yes
else
    ac_enable_swan=no
fi
])

dnl Add option for --enable-mpi
ac_have_mpi=no
AC_ARG_ENABLE(mpi, AS_HELP_STRING([--enable-mpi],[Build with parallel support]),[
if test "x$enable_mpi" = "xyes" ; then
    ac_build_mpi=yes
else
    ac_build_mpi=no
fi
])

dnl If we are using MPI, then check for the compiler wrapper
if test "x$ac_build_mpi" = "xyes" ; then
    AC_LANG_PUSH([Fortran])
    AX_MPI([ac_have_mpi=yes],[AC_MSG_ERROR([Could not find MPI library and compilers])])
    AC_LANG_POP()
fi
  
dnl Add option for --enable-netcdf
AC_ARG_ENABLE(netcdf, AS_HELP_STRING([--enable-netcdf], [Build with netCDF support]))

dnl If we are using netcdf, check the configuration here
nc_has_fortran=0
if test "x$enable_netcdf" = "xyes" ; then
  AC_LANG_PUSH([C])
  AC_CHECK_LIB([netcdf],[nc_open])
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
      AC_MSG_ERROR([netcdf-fortran is required to enable netcdf. Try configuring with --disable-netcdf])
  fi

  AC_LANG_POP()
fi

dnl If building swan, we need perl
AS_IF([test x$ac_enable_swan = "xyes"],
    [
     AC_CHECK_PROG([HAS_PERL],[perl],[yes],[no])
     AS_IF([test x$HAS_PERL != xyes],AC_MSG_ERROR([Perl required to build SWAN source code]),[])
    ],[])

dnl To generate the version file, we need python
AC_CHECK_PROG([HAS_PYTHON],[python3],[yes],[no])
AS_IF([test x$HAS_PYTHON != xyes],AC_MSG_ERROR([Python require to generate the version information]),[])

dnl Set the conditional values in all Makefile.am files 
AM_CONDITIONAL([BUILD_MPI], [test x$ac_have_mpi = xyes])
AM_CONDITIONAL([BUILD_NETCDF], [test x$nc_has_fortran = xyes])
AM_CONDITIONAL([BUILD_SWAN], [test x$ac_enable_swan = xyes])

dnl Generate the SWAN source files via the preprocessor
AS_IF([test x$ac_enable_swan = "xyes"],
    [
     AC_CONFIG_COMMANDS([adcswan], 
        [
         perl thirdparty/swan/switch.pl -unix -adcirc thirdparty/swan/*.ftn thirdparty/swan/*.ftn90
         mkdir -p work/lib/libswan_serial_source
         mv thirdparty/swan/*.f thirdparty/swan/*.f90 work/lib/libswan_serial_source/.
        ]
     )
    ],[]
)
AS_IF([test x$ac_enable_swan = "xyes" && test x$ac_have_mpi = xyes],
    [
     AC_CONFIG_COMMANDS([padcswan], 
        [
         perl thirdparty/swan/switch.pl -unix -adcirc -pun thirdparty/swan/*.ftn thirdparty/swan/*.ftn90
         mkdir -p work/lib/libswan_parallel_source
         mv thirdparty/swan/*.f thirdparty/swan/*.f90 work/lib/libswan_parallel_source/.
        ]
     )
    ],[]
)

dnl Create the adcirc version source file
AC_CONFIG_COMMANDS([version],
  [scripts/adcirc_version.py --create-version-file --directory $ac_top_srcdir]
)

dnl Set the Makefiles to be created
AC_CONFIG_FILES([
  Makefile
  work/Makefile
  work/lib/Makefile
  work/lib/metis/Makefile
  work/lib/libswan_serial1/Makefile
  work/lib/libadcswan/Makefile
  work/lib/libswan_serial2/Makefile
  work/lib/libswan_parallel1/Makefile
  work/lib/libpadcswan/Makefile
  work/lib/libswan_parallel2/Makefile
  work/bin/Makefile
  work/bin/adcirc/Makefile
  work/bin/adcswan/Makefile
  work/bin/adcprep/Makefile
  work/bin/padcirc/Makefile
  work/bin/padcswan/Makefile
  work/bin/aswip/Makefile
])

dnl Write the output
AC_OUTPUT
