F77 = $(MPIFC)
FC = $(MPIFC)

VPATH = $(top_srcdir) \
		$(top_srcdir)/src \
		$(top_srcdir)/prep \
		$(top_srcdir)/wind \
		$(top_srcdir)/thirdparty/KDTREE2

bin_PROGRAMS = padcirc

padcirc_SOURCES = \
	mkdir.c \
	version.F \
	sizes.F \
	constants.F \
	kdtree2.F \
	global.F \
	boundaries.F \
	global_3dvs.F \
	messenger.F \
	netcdf_error.F90 \
	hashtable.F \
	mesh.F \
	subgridLookup.F \
	internaltide.F \
	nodalattr.F \
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
	writer.F \
	write_output.F \
	wetdry.F \
	rs2.F \
	adcirc.F \
	weir_boundary.F \
	gl2loc_mapping.F \
	couple2baroclinic3D.F \
	transport.F \
	vsmy.F \
	timestep.F \
	cstart.F \
	hstart.F \
	read_input.F \
	control.F \
    driver.F 		

MPI_FLAGS=-DCMPI
ARCH_FALGS=-DLINUX -DCSCA
if HAVE_NETCDF
NETCDF_FLAGS = -DADCNETCDF -DHAVE_NETCDF4 -DNETCDF_CAN_DEFLATE
AM_FFLAGS = $(NETCDF_FLAGS) $(FCFLAGS) $(MPI_FLAGS) $(ARCH_FLAGS)
AM_FCFLAGS = $(NETCDF_FLAGS) $(MPI_FLAGS) $(ARCH_FLAGS)
padcirc_LDADD = -lnetcdf -lnetcdff 
else
AM_FFLAGS = $(FCFLAGS) $(MPI_FLAGS) $(ARCH_FLAGS)
AM_FCFLAGS = $(ARCH_FLAGS)
endif

CLEANFILES = *.mod
