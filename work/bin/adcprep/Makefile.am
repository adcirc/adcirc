F77 = $(FC)

VPATH = $(top_srcdir) \
		$(top_srcdir)/src \
		$(top_srcdir)/prep \
		$(top_srcdir)/wind \
		$(top_srcdir)/thirdparty/KDTREE2

if BUILD_MPI

bin_PROGRAMS = adcprep

adcprep_SOURCES = \
    mkdir.c version.F sizes.F constants.F kdtree2.F global.F boundaries.F hashtable.F \
    mesh.F netcdf_error.F90 subgridLookup.F internaltide.F nodalattr.F vew1d.F \
	  global_3dvs.F vortex.F owiwind.F rs2.F owi_ice.F owiwind_netcdf.F nws08.F90 wind.F \
    presizes.F pre_global.F metis.F subprep.F adcprep.F decomp.F prep_weir.F \
    itpackv.F harm.F read_global.F subdomain.F quadrature.F sponge_layer.F \
	  momentum.F gwce.F wetdry.F netcdfio.F prep.F interp.F machdep.F \
	  couple2baroclinic3D.F gl2loc_mapping.F 

adcprep_LDADD = $(top_builddir)/work/lib/metis/libmetis.a

ARCH_FLAGS = -DLINUX -DCSCA

if BUILD_NETCDF
NETCDF_FLAGS = -DADCNETCDF -DHAVE_NETCDF4 -DNETCDF_CAN_DEFLATE
AM_FFLAGS = $(NETCDF_FLAGS) $(FCFLAGS) $(ARCH_FLAGS)
AM_FCFLAGS = $(NETCDF_FLAGS)
adcprep_LDADD += -lnetcdf -lnetcdff 
else
AM_FFLAGS = $(FCFLAGS) $(ARCH_FLAGS)
AM_FCFLAGS = 
endif

if BUILD_SWAN
AM_FFLAGS += -DADCSWAN
AM_FCFLAGS += -DADCSWAN
endif

endif

CLEANFILES = *.mod
