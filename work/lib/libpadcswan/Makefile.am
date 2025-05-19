include ../../flags.am

F77 = $(MPIFC)
FC = $(MPIFC)

VPATH = $(top_srcdir) \
		$(top_srcdir)/src \
		$(top_srcdir)/wind \
		$(top_srcdir)/prep \
		$(top_srcdir)/thirdparty/KDTREE2

internal_includes = -I$(top_srcdir)/work/lib/libswan_parallel1 
AM_FFLAGS += $(internal_includes) -DCSWAN -DCMPI
AM_FCFLAGS += $(internal_includes) -DCSWAN -DCMPI
ARFLAGS = cr

if BUILD_SWAN

if BUILD_MPI

noinst_LIBRARIES = libpadcswan.a

libpadcswan_a_SOURCES = \
	  mkdir.c version.F sizes.F constants.F kdtree2.F global.F boundaries.F \
    hashtable.F global_3dvs.F messenger.F mesh.F netcdf_error.F90 subgridLookup.F \
	  internaltide.F nodalattr.F vew1d.F harm.F vortex.F owiwind.F \
	  owiwind_netcdf.F nws08.F90 wind.F owiwind.F rs2.F owi_ice.F itpackv.F \
	  globalio.F subdomain.F quadrature.F sponge_layer.F momentum.F \
	  gwce.F wetdry.F control.F netcdfio.F writer.F write_output.F couple2swan.F \
	  gl2loc_mapping.F couple2baroclinic3D.F

endif

endif

CLEANFILES = *.mod
