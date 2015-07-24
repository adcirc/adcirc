NETCDFPATH  := /usr
XDMFPATH    := /home/jason/projects/XDMF/Code/latestCode
XDMFLIBPATH := /home/jason/projects/XDMF/Code/testLatest
HDF5PATH    := /home/jason/projects/XDMF/Code/hdf5-1.8.11/hdf5
ASGSPATH    := /home/jason/asgs/trunk
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
#
OBJ := adcmesh.o control.o nodalattr.o asgsio.o
MODS := control.mod adcmesh.mod nodalattr.mod asgsio.mod
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
control.mod : control.f90 adcmesh.mod
	$(FC) $(FFLAGS) $(INCLUDES) $(LIBS) -c control.f90 $(LDFLAGS)
#
nodalattr.mod : $(ASGSPATH)/util/nodalattr/nodalattr.f90
	$(FC) $(FFLAGS) $(INCLUDES) $(LIBS) -c $(ASGSPATH)/util/nodalattr/nodalattr.f90 $(LDFLAGS)
