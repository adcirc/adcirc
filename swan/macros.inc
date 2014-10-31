##############################################################################
# IA32_Intel/x86-64_Intel:	Intel Pentium with Linux using Intel compiler 11.
##############################################################################
F90_SER = gfortran
F90_OMP = gfortran
F90_MPI = mpif90
FLAGS_OPT = -O2
FLAGS_MSC = -W0 -assume byterecl -traceback -diag-disable 8291 -diag-disable 8293
FLAGS90_MSC = $(FLAGS_MSC)
FLAGS_SER =
FLAGS_OMP = -openmp
FLAGS_MPI =
INCS_SER =
INCS_OMP =
INCS_MPI =
LIBS_SER =
LIBS_OMP =
LIBS_MPI =
O_DIR = ../work/odir4/
OUT = -o 
EXTO = o
MAKE = make
RM = rm -f
swch = -unix
