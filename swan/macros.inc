##############################################################################
# IA32_Intel/x86-64_Intel:	Intel Pentium with Linux using Intel compiler 11.
##############################################################################
F90_SER = ifort
F90_OMP = ifort
F90_MPI = mpif90
FLAGS_OPT = -O2
FLAGS_MSC = -W0 -assume byterecl # -traceback -diag-disable 8290 -diag-disable 8291 -diag-disable 8293
FLAGS90_MSC = $(FLAGS_MSC)
FLAGS_SER =
FLAGS_OMP = -openmp
FLAGS_MPI =
NETCDFROOT =
ifneq ($(NETCDFROOT),)
  INCS_SER = -I$(NETCDFROOT)/include
  INCS_OMP = -I$(NETCDFROOT)/include
  INCS_MPI = -I$(NETCDFROOT)/include
  LIBS_SER = -L$(NETCDFROOT)/lib -lnetcdf -lnetcdff
  LIBS_OMP = -L$(NETCDFROOT)/lib -lnetcdf -lnetcdff
  LIBS_MPI = -L$(NETCDFROOT)/lib -lnetcdf -lnetcdff
  NCF_OBJS = nctablemd.o agioncmd.o swn_outnc.o
else
  INCS_SER =
  INCS_OMP =
  INCS_MPI =
  LIBS_SER =
  LIBS_OMP =
  LIBS_MPI =
  NCF_OBJS =
endif
O_DIR = ../work/odir4/
OUT = -o 
EXTO = o
MAKE = make
RM = rm -f
ifneq ($(NETCDFROOT),)
  swch = -unix -impi -netcdf
else
  swch = -unix -impi
endif
