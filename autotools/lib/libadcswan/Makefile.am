include ../../flags.am

F77 = $(FC)

VPATH = $(top_srcdir) \
		$(top_srcdir)/src \
		$(top_srcdir)/wind \
		$(top_srcdir)/thirdparty/KDTREE2

noinst_LIBRARIES = libadcswan.a

libadcswan_a_SOURCES = \
	  version.F sizes.F constants.F kdtree2.F global.F boundaries.F \
      hashtable.F mesh.F netcdf_error.F90 subgridLookup.F internaltide.F \
	  nodalattr.F vew1d.F global_3dvs.F harm.F vortex.F owiwind.F \
	  owiwind_netcdf.F wind.F owiwind.F rs2.F owi_ice.F itpackv.F \
	  globalio.F subdomain.F quadrature.F sponge_layer.F momentum.F \
	  gwce.F wetdry.F control.F netcdfio.F write_output.F couple2swan.F \
	  couple2baroclinic3D.F gl2loc_mapping.F couple2baroclinic3D.F

internal_includes = -I$(top_srcdir)/autotools/lib/libswan_serial1 

AM_FFLAGS += $(internal_includes) -DCSWAN
AM_FCFLAGS += $(internal_includes) -DCSWAN

CLEANFILES = *.mod
