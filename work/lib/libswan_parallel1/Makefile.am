VPATH = $(top_srcdir)/work/lib/libswan_parallel_source

F77 = $(MPIFC)
FC = $(MPIFC)

ARFLAGS = cr

if BUILD_SWAN

if BUILD_MPI

noinst_LIBRARIES = libswan_parallel1.a

libswan_parallel1_a_SOURCES = \
    swmod1.f swmod2.f SwanSpectPart.f m_constants.f90 m_fileio.f90 \
    serv_xnl4v5.f90 mod_xnl4v5.f90 SwanGriddata.f90 SwanGridobjects.f90 \
    SwanCompdata.f90 couple2adcirc.f90 SdsBabanin.f90 SwanIEM.f90 \
	  SwanBraggScat.f90 SwanQCM.f90

endif

endif
