##############################################################################
# IA32_GNU:		Intel Pentium with Linux using GNU compiler gfortran.
##############################################################################
F90_SER = gfortran
F90_OMP = gfortran
F90_MPI = mpif90
FLAGS_OPT = -O
FLAGS_MSC = -w -fno-second-underscore
FLAGS90_MSC = $(FLAGS_MSC) -ffree-line-length-none
FLAGS_DYN =
FLAGS_SER =
FLAGS_OMP = -fopenmp
FLAGS_MPI =
NETCDFROOT =
ifneq ($(NETCDFROOT),)
  INCS_SER = -I$(NETCDFROOT)/include
  INCS_OMP = -I$(NETCDFROOT)/include
  INCS_MPI = -I$(NETCDFROOT)/include
  LIBS_SER = -L$(NETCDFROOT)/lib -lnetcdf -lnetcdff
  LIBS_OMP = -L$(NETCDFROOT)/lib -lnetcdf -lnetcdff -static-libgcc
  LIBS_MPI = -L$(NETCDFROOT)/lib -lnetcdf -lnetcdff
  NCF_OBJS = nctablemd.o agioncmd.o swn_outnc.o
else
  INCS_SER =
  INCS_OMP =
  INCS_MPI =
  LIBS_SER =
  LIBS_OMP = -static-libgcc
  LIBS_MPI =
  NCF_OBJS =
endif
O_DIR = ../work/odir4/
OUT = -o 
EXTO = o
MAKE = make
RM = rm -f
ifneq ($(NETCDFROOT),)
  swch = -unix -netcdf
else
  swch = -unix
endif
