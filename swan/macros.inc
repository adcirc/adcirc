##############################################################################
# IA32_Intel/x86-64_Intel:	Intel Pentium with Linux using Intel compiler 9.1.
##############################################################################
F90_SER = ifort
F90_OMP = ifort
F90_MPI = mpif90
FLAGS_OPT = -O2 -ipo -xW -mp1
FLAGS_MSC = -W0 -assume byterecl -traceback
FLAGS_F90 =
FLAGS_SER =
FLAGS_OMP = -openmp -assume cc_omp -fpp2
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
swch = -unix -f95 -timg -impi
