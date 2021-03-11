SHELL=/bin/sh
# If it is invoked by the command line
# #       make -f makefile
#
SRCS=estofs_tide_fac.f

OBJS=estofs_tide_fac.o

# Tunable parameters
# FC		Name of the fortran compiling system to use
# LDFLAGS	Flags to the loader
# LIBS		List of libraries
# CMD		Name of the executable
# PROFLIB	Library needed for profiling

FC =            ifort
#use for extensive compilation error reporting
LDFLAGS =	-mkl  # ifort mkl library

CMD =           estofs_tide_fac
PROFLIB =       -lprof

# To perform the default compilation, use the first line
# To compile with flowtracing turned on, use the second line
# To compile giving profile additonal information, use the third line
# WARNING:  SIMULTANEOUSLY PROFILING AND FLOWTRACING IS NOT RECOMMENDED
FFLAGS =        -O3 -i4 $(INCS)

# Lines from here on down should not need to be changed.  They are the
# actual rules which make uses to build a.out.

all:		$(CMD)

$(CMD):		$(OBJS)
	$(FC) $(LDFLAGS)  -o $(@) $(OBJS) $(LIBS)
	rm -f $(OBJS)
#	cp -p $(CMD) ../../exec/
# Make the profiled version of the command and call it a.out.prof

#$(FC) $(LDFLAGS)  -o $(@) $(OBJS) $(PROFLIB) $(LIBS)

clean:
	rm -f $(OBJS)

clobber:	clean
	rm -f $(CMD) $(CMD).prof

void:	clobber
	rm -f $(SRCS) makefile

