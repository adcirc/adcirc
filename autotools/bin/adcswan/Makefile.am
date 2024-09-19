include ../../flags.am

F77 = $(FC)

VPATH = $(top_srcdir)/src

if BUILD_SWAN

bin_PROGRAMS = adcswan

adcswan_SOURCES = \
      adcirc.F weir_boundary.F read_input.F cstart.F hstart.F timestep.F \
      vsmy.F transport.F driver.F

adcswan_LDADD = \
				$(top_builddir)/autotools/lib/libadcswan/libadcswan.a \
				$(top_builddir)/autotools/lib/libswan_serial2/libswan_serial2.a \
				$(top_builddir)/autotools/lib/libswan_serial1/libswan_serial1.a 

internal_includedirs = -I../../lib/libswan_serial1 -I../../lib/libswan_serial2 -I../../lib/libadcswan

AM_FFLAGS += $(internal_includedirs) -DCSWAN
AM_FCFLAGS += $(internal_includedirs) -DCSWAN

if BUILD_NETCDF
adcswan_LDADD += -lnetcdf -lnetcdff
endif

endif

CLEANFILES = *.mod
