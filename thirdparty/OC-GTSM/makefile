SHELL:=/bin/bash

FC = mpif90 -cpp
O_DIR =  
LIBHOME    = ../
HDF5HOME   = $(HDF5_HOME)
NETCDFHOME = $(NETCDF_HOME)
GSWHOME       := $(LIBHOME)GSW-Fortran/build/
DATETIMEHOME  := $(LIBHOME)datetime-fortran/build/
incdir = -I$(DATETIMEHOME)include -I$(NETCDFHOME)/include -I$(GSWHOME)gsw/
libdir = -L$(NETCDFHOME)/lib -L$(HDF5HOME) -L$(GSWHOME) -L$(DATETIMEHOME)lib/
FFLAGS = -O2 -c $(incdir) -fbacktrace #-g -check bounds 
LFLAGS = -O2 $(libdir) -lnetcdf -lnetcdff -lgsw -ldatetime
LINK = $(FC)
TARGET = OGCM_DL.a

$(O_DIR):
	mkdir -p $@

SRC  =  OGCM_DL.f90
	
OBJ:=   $(patsubst %.f90, $(O_DIR)%.o, $(SRC) )

all: $(TARGET)

clean: 
	-rm -rf $(O_DIR)*.o

$(TARGET): $(OBJ)
	$(LINK) -o $@ $(OBJ) $(LFLAGS) 

$(O_DIR)%.o : %.f90
	$(FC) $(FFLAGS) $< -o $@
