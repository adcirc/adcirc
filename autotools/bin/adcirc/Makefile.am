include ../../flags.am

F77 = $(FC)

VPATH = $(top_srcdir) \
		$(top_srcdir)/src \
		$(top_srcdir)/wind \
		$(top_srcdir)/thirdparty/KDTREE2

bin_PROGRAMS = adcirc

adcirc_SOURCES = \
	version.F sizes.F constants.F kdtree2.F global.F netcdf_error.F90 hashtable.F \
	boundaries.F mesh.F subgridLookup.F internaltide.F nodalattr.F global_3dvs.F \
	itpackv.F quadrature.F sponge_layer.F subdomain.F owi_ice.F owiwind.F \
	owiwind_netcdf.F vortex.F wind.F vew1d.F momentum.F gwce.F harm.F netcdfio.F \
	globalio.F write_output.F wetdry.F rs2.F adcirc.F weir_boundary.F couple2baroclinic3D.F \
	transport.F vsmy.F timestep.F cstart.F hstart.F read_input.F control.F driver.F 		

if BUILD_NETCDF
adcirc_LDADD = -lnetcdf -lnetcdff 
endif

CLEANFILES = *.mod
