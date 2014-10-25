NETCDFPATH  := /usr
XDMFPATH    := /home/jason/projects/XDMF/Code/latestCode
XDMFLIBPATH := /home/jason/projects/XDMF/Code/testLatest
ASGSPATH    := /home/jason/asgs/trunk
LDFLAGS     := -lXdmfCore -lXdmfUtils -lXdmf -lnetcdff
LIBS        := -L$(XDMFLIBPATH) -L$(NETCDFPATH)/lib
INCLUDES    := -I$(XDMFPATH)/utils -I$(ASGSPATH)/output -I$(NETCDFPATH)/include
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
OBJ := adcmesh.o control.o nodalattr.o
MODS := control.mod adcmesh.mod nodalattr.mod
#
# targets
all : adcirc2xdmf xdmf2adcirc
#
clean :
	rm -f *.o *.mod *.x
#
adcirc2xdmf : $(MODS)
	$(FC) $(FFLAGS) $(INCLUDES) $(LIBS) -o $@.x $@.f90 $(OBJ) $(LDFLAGS)
#
xdmf2adcirc : $(MODS)
	$(FC) $(FFLAGS) $(INCLUDES) $(LIBS) -o $@.x $@.f90 $(OBJ) $(LDFLAGS)
#
adcmesh.mod : $(ASGSPATH)/output/adcmesh.f90
	$(FC) $(FFLAGS) $(INCLUDES) $(LIBS) -c $(ASGSPATH)/output/adcmesh.f90 $(LDFLAGS)
#
control.mod : control.f90 adcmesh.mod
	$(FC) $(FFLAGS) $(INCLUDES) $(LIBS) -c control.f90 $(LDFLAGS)
#
nodalattr.mod : $(ASGSPATH)/util/nodalattr/nodalattr.f90
	$(FC) $(FFLAGS) $(INCLUDES) $(LIBS) -c $(ASGSPATH)/util/nodalattr/nodalattr.f90 $(LDFLAGS)
