include ../../flags.am

F77 = $(FC)

VPATH = $(top_srcdir) \
		$(top_srcdir)/src \
		$(top_srcdir)/wind \
		$(top_srcdir)/thirdparty/KDTREE2

bin_PROGRAMS = aswip

aswip_SOURCES = \
      sizes.F \
      constants.F \
      kdtree2.F \
      global.F \
      hashtable.F \
      boundaries.F \
      mesh.F \
			netcdf_error.F90 \
      subgridLookup.F \
      internaltide.F \
      nodalattr.F \
      global_3dvs.F \
      gl2loc_mapping.F \
      owiwind.F \
			owiwind_netcdf.F \
      owi_ice.F \
      vortex.F \
      nws08.F90 \
      wind.F \
      aswip.F

if BUILD_NETCDF
aswip_LDADD = -lnetcdf -lnetcdff 
endif

CLEANFILES = *.mod
