Chris Massey, USACE-ERDC-CHL, Vicksburg, MS 39180
May 6, 2011

Made a couple of modifications to the packages called by metis to get
it to work on a Windows based machine.  Added the gettimeofday.h module
which is not available on a Windows machine.  A person wanting to compile
metis for a Windows based PC would select all the routines in the Lib
folder and would also include the gettimeofday.h routine in the Windows
folder and would use a -WINDOWS pre-processor flag in addition to the
other pre-processor flags in the normal makefile.  NOTE:  Only use the
gettimeofday.h for Windows machines, Linux/Unix machines typically already
have it installed.

These changes to metis do not change the usual building or functionality
of metis on Linux/Unix machines.


Pre-processor definitions to use on a Window's PC for compiling metis:
NDEBUG;_LIB;__VC__;WINDOWS
