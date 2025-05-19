include ../../flags.am

F77 = $(MPIFC)
FC = $(MPIFC)

VPATH = $(top_srcdir) \
		$(top_srcdir)/src \
		$(top_srcdir)/prep \
		$(top_srcdir)/wind \
		$(top_srcdir)/thirdparty/KDTREE2

if BUILD_MPI

bin_PROGRAMS = padcirc

padcirc_SOURCES = \
	mkdir.c version.F sizes.F constants.F kdtree2.F global.F boundaries.F \
	global_3dvs.F messenger.F netcdf_error.F90 hashtable.F mesh.F subgridLookup.F \
	internaltide.F nodalattr.F itpackv.F quadrature.F sponge_layer.F subdomain.F \
	owi_ice.F owiwind.F owiwind_netcdf.F vortex.F nws08.F90 wind.F vew1d.F momentum.F \
	gwce.F harm.F netcdfio.F globalio.F writer.F write_output.F wetdry.F \
	rs2.F adcirc.F weir_boundary.F gl2loc_mapping.F couple2baroclinic3D.F \
	astronomic.F90 sun.F90 moon.F90 sun_moon_system.F90 ephemerides.F90 tidalpotential.F90 \
	transport.F vsmy.F timestep.F cstart.F hstart.F read_input.F control.F \
  driver.F 		

AM_FCFLAGS += -DCMPI
AM_FFLAGS += -DCMPI
if BUILD_NETCDF
padcirc_LDADD = -lnetcdf -lnetcdff 
endif

endif

CLEANFILES = *.mod
