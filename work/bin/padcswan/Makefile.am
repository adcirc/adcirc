include ../../flags.am

F77 = $(MPIFC)
FC = $(MPIFC)

VPATH = $(top_srcdir)/src

if BUILD_MPI 

if BUILD_SWAN

bin_PROGRAMS = padcswan

padcswan_SOURCES = \
      adcirc.F weir_boundary.F astronomic.F90 sun.F90 moon.F90 \
      sun_moon_system.F90 ephemerides.F90 tidalpotential.F90 \
      timestep.F cstart.F hstart.F read_input.F vsmy.F transport.F \
      driver.F

padcswan_LDADD = \
				$(top_builddir)/work/lib/libpadcswan/libpadcswan.a \
				$(top_builddir)/work/lib/libswan_parallel2/libswan_parallel2.a \
				$(top_builddir)/work/lib/libswan_parallel1/libswan_parallel1.a 

internal_includedirs = -I../../lib/libswan_parallel1 -I../../lib/libswan_parallel2 -I../../lib/libpadcswan

AM_FFLAGS += $(internal_includedirs) -DCSWAN -DCMPI
AM_FCFLAGS += $(internal_includedirs) -DCSWAN -DCMPI

if BUILD_NETCDF
padcswan_LDADD += -lnetcdf -lnetcdff
endif

endif

endif

CLEANFILES = *.mod
