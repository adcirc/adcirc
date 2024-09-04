F77 = $(FC)

VPATH = $(top_srcdir) \
		$(top_srcdir)/src \
		$(top_srcdir)/wind \
		$(top_srcdir)/thirdparty/KDTREE2

bin_PROGRAMS = adcirc

adcirc_SOURCES = \
	version.F \
	sizes.F \
	constants.F \
	kdtree2.F \
	global.F \
	netcdf_error.F90 \
	hashtable.F \
	boundaries.F \
	mesh.F \
	subgridLookup.F \
	internaltide.F \
	nodalattr.F \
	global_3dvs.F \
	itpackv.F \
	quadrature.F \
	sponge_layer.F \
	subdomain.F \
	owi_ice.F \
	owiwind.F \
	owiwind_netcdf.F \
	vortex.F \
	wind.F \
	vew1d.F \
	momentum.F \
	gwce.F \
	harm.F \
	netcdfio.F \
	globalio.F \
	write_output.F \
	wetdry.F \
	rs2.F \
	adcirc.F \
	weir_boundary.F \
	couple2baroclinic3D.F \
	transport.F \
	vsmy.F \
	timestep.F \
	cstart.F \
	hstart.F \
	read_input.F \
	control.F \
    driver.F 		

ARCH_FLAGS=-DLINUX -DCSCA
if HAVE_NETCDF
NETCDF_FLAGS = -DADCNETCDF -DHAVE_NETCDF4 -DNETCDF_CAN_DEFLATE
AM_FFLAGS = $(NETCDF_FLAGS) $(FCFLAGS) $(ARCH_FLAGS)
AM_FCFLAGS = $(NETCDF_FLAGS) $(ARCH_FLAGS)
adcirc_LDADD = -lnetcdf -lnetcdff 
else
AM_FFLAGS = $(FCFLAGS) $(ARCH_FLAGS)
AM_FCFLAGS = $(ARCH_FLAGS)
endif

CLEANFILES = *.mod
