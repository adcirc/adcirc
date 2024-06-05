################################################################################
# IA32_Intel/x86-64_Intel:	Intel Core/Xeon with Linux using Intel compiler.
################################################################################
F90_SER = ifort
F90_OMP = ifort
# for older compiler versions, use mpif90
F90_MPI = mpiifort
FLAGS_OPT = -O2
FLAGS_MSC = -W0 -assume byterecl -traceback -diag-disable 8290 -diag-disable 8291 -diag-disable 8293
FLAGS90_MSC = $(FLAGS_MSC)
FLAGS_DYN = -fPIC
FLAGS_SER =
FLAGS_OMP = -qopenmp
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
O_DIR = ../../work/odir4/
BOU_OBJ = $(O_DIR)boundaries.o
ifneq ("$(wildcard $(BOU_OBJ))","")
  MSG_OBJS = $(O_DIR)mkdir.o $(O_DIR)boundaries.o $(O_DIR)sizes.o $(O_DIR)global.o $(O_DIR)global_3dvs.o $(O_DIR)version.o $(O_DIR)messenger.o
else
  MSG_OBJS = $(O_DIR)mkdir.o $(O_DIR)sizes.o $(O_DIR)global.o $(O_DIR)global_3dvs.o $(O_DIR)version.o $(O_DIR)messenger.o
endif
OUT = -o 
EXTO = o
MAKE = make
RM = rm -f
ifneq ($(NETCDFROOT),)
  swch = -unix -impi -netcdf
else
  swch = -unix -impi
endif
