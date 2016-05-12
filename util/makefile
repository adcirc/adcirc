NETCDFPATH  := /usr
XDMFPATH    := /home/jason/projects/XDMF/Code/latestCode
XDMFLIBPATH := /home/jason/projects/XDMF/Code/testLatest
HDF5PATH    := /home/jason/projects/XDMF/Code/hdf5-1.8.11/hdf5
ASGSPATH    := /home/jason/asgs/master
LIBS        := -L$(XDMFLIBPATH) -L$(NETCDFPATH)/lib -L$(HDF5PATH)/lib
INCLUDES    := -I$(XDMFPATH)/utils -I$(ASGSPATH)/output -I$(NETCDFPATH)/include
LDFLAGS     := -lXdmfCore -lXdmfUtils -lXdmf -lnetcdff -Wl,-rpath=$(HDF5PATH)/lib
# 
# specify compiler=gfortran on the make command line
ifeq ($(compiler),gfortran)
   FC := gfortran
   FFLAGS := -cpp -ffree-line-length-none
   ifeq ($(DEBUG),full)
      FFLAGS := -cpp -g -O0 -fbacktrace -fbounds-check -ffpe-trap=zero,invalid,underflow,overflow,denormal #-Wall
   endif
endif
#
# If netcdf libraries are available on the target platform, 
# use NETCDF=enable on the make command line; this is not needed
# for xdmf2adcirc or adcirc2xdmf, but it is needed by their dependencies.
ifeq ($(NETCDF),enable)
   FFLAGS      := $(FFLAGS) -DASGSNETCDF
   LDFLAGS     := $(LDFLAGS) -lnetcdff
   ifeq ($(MACHINENAME),jason-desktop)
      NETCDFPATH  := /usr
   endif
   ifeq ($(MACHINENAME),hatteras)
      NETCDFPATH  := /usr/share/Modules/software/RHEL-6.5/netcdf/netcdf-4.1.3_intel-14.0.3
   endif
   # jgf20150616: Added netcdf paths for stampede at TACC.
   ifeq ($(MACHINENAME),stampede)
      NETCDFPATH := /opt/apps/intel13/netcdf/4.3.2/x86_64
   endif
   # jgf20150720: Added compiler flags for queenbee at LONI.
   ifeq ($(MACHINENAME),queenbee)
      NETCDFPATH  := /usr/local/packages/netcdf/4.2.1.1/INTEL-140-MVAPICH2-2.0
   endif
   # jgf20150817: Added compiler flags for garnet at ERDC.
   ifeq ($(MACHINENAME),garnet)
      NETCDFPATH  := /opt/cray/netcdf/4.3.0/pgi/121
   endif
   FFLAGS      := $(FFLAGS) -DASGSNETCDF
   INCLUDES    := $(INCLUDES) -I. -I$(NETCDFPATH)/include
   LIBS        := $(LIBS) -L$(NETCDFPATH)/lib
endif
ifeq ($(NETCDF4),enable)
   FFLAGS      := $(FFLAGS) -DHAVE_NETCDF4
endif
ifeq ($(NETCDF4_COMPRESSION),enable)
   FFLAGS      := $(FFLAGS) -DNETCDF_CAN_DEFLATE
endif
#
#
OBJ := adcmesh.o control.o nodalattr.o asgsio.o adcircdata.o
MODS := control.mod adcmesh.mod nodalattr.mod asgsio.mod adcircdata.mod
#
# targets
all : adcirc2xdmf xdmf2adcirc
#
clean :
	rm -f *.o *.mod *.x $(ASGSPATH)/output/adcmesh.o $(ASGSPATH)/output/asgsio.o $(ASGSPATH)/util/nodalattr/nodalattr.o $(ASGSPATH)/output/adcmesh.mod $(ASGSPATH)/output/asgsio.mod $(ASGSPATH)/util/nodalattr/nodalattr.mod  
#
adcirc2xdmf : $(MODS)
	$(FC) $(FFLAGS) $(INCLUDES) $(LIBS) -o $@.x $@.f90 $(OBJ) $(LDFLAGS)
#
xdmf2adcirc : $(MODS)
	$(FC) $(FFLAGS) $(INCLUDES) $(LIBS) -o $@.x $@.f90 $(OBJ) $(LDFLAGS)
#
adcmesh.mod : $(ASGSPATH)/output/adcmesh.f90 asgsio.mod
	$(FC) $(FFLAGS) $(INCLUDES) $(LIBS) -c $(ASGSPATH)/output/adcmesh.f90 $(LDFLAGS)
#
asgsio.mod : $(ASGSPATH)/output/asgsio.f90
	$(FC) $(FFLAGS) $(INCLUDES) $(LIBS) -c $(ASGSPATH)/output/asgsio.f90 $(LDFLAGS)
#
adcircdata.mod : $(ASGSPATH)/output/adcircdata.f90
	$(FC) $(FFLAGS) $(INCLUDES) $(LIBS) -c $(ASGSPATH)/output/adcircdata.f90 $(LDFLAGS)
#
control.mod : control.f90 adcmesh.mod
	$(FC) $(FFLAGS) $(INCLUDES) $(LIBS) -c control.f90 $(LDFLAGS)
#
nodalattr.mod : $(ASGSPATH)/util/nodalattr/nodalattr.f90
	$(FC) $(FFLAGS) $(INCLUDES) $(LIBS) -c $(ASGSPATH)/util/nodalattr/nodalattr.f90 $(LDFLAGS)
