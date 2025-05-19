include ../../flags.am

F77 = $(FC)

VPATH = $(top_srcdir)/src

if BUILD_SWAN

bin_PROGRAMS = adcswan

adcswan_SOURCES = \
      adcirc.F weir_boundary.F astronomic.F90 sun.F90 moon.F90 \
      sun_moon_system.F90 ephemerides.F90 tidalpotential.F90 \
      timestep.F cstart.F hstart.F read_input.F vsmy.F transport.F \
      driver.F

adcswan_LDADD = \
				$(top_builddir)/work/lib/libadcswan/libadcswan.a \
				$(top_builddir)/work/lib/libswan_serial2/libswan_serial2.a \
				$(top_builddir)/work/lib/libswan_serial1/libswan_serial1.a 

internal_includedirs = -I../../lib/libswan_serial1 -I../../lib/libswan_serial2 -I../../lib/libadcswan

AM_FFLAGS += $(internal_includedirs) -DCSWAN
AM_FCFLAGS += $(internal_includedirs) -DCSWAN

if BUILD_NETCDF
adcswan_LDADD += -lnetcdf -lnetcdff
endif

endif

CLEANFILES = *.mod
