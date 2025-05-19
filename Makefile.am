ACLOCAL_AMFLAGS = -I m4
SUBDIRS = work

adcirc:
	cd work/bin/adcirc && $(MAKE) $(AM_MAKEFLAGS) adcirc 

aswip:
	cd work/bin/aswip && $(MAKE) $(AM_MAKEFLAGS) aswip

if BUILD_MPI

adcprep:
	cd work/bin/adcprep && $(MAKE) $(AM_MAKEFLAGS) adcprep

padcirc:
	cd work/bin/padcirc && $(MAKE) $(AM_MAKEFLAGS) padcirc

endif

if BUILD_SWAN

adcswan:
	cd work/bin/adcswan && $(MAKE) $(AM_MAKEFLAGS) adcswan
	
if BUILD_MPI

padcswan:
	cd work/bin/padcswan && $(MAKE) $(AM_MAKEFLAGS) padcswan

endif

endif


.PHONY: adcirc adcprep padcirc adcswan padcswan aswip 
