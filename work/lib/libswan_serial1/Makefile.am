VPATH = $(top_srcdir)/work/lib/libswan_serial_source

ARFLAGS = cr

if BUILD_SWAN

noinst_LIBRARIES = libswan_serial1.a

libswan_serial1_a_SOURCES = \
    swmod1.f swmod2.f SwanSpectPart.f m_constants.f90 m_fileio.f90 \
    serv_xnl4v5.f90 mod_xnl4v5.f90 SwanGriddata.f90 SwanGridobjects.f90 \
    SwanCompdata.f90 couple2adcirc.f90 SdsBabanin.f90 SwanIEM.f90 \
	  SwanBraggScat.f90 SwanQCM.f90

endif
