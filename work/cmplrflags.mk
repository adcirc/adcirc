# SRCDIR is set in makefile or on the compile line
INCDIRS := -I . -I $(SRCDIR)/prep

INCDIRS += -I$(NETCDFHOME)/include/  
LIBDIRS += -L$(NETCDFHOME)/lib/

INCDIRS += -I$(HDF5HOME)/include/   
LIBDIRS += -L$(HDF5HOME)/lib/

LIBDIRS += -DALL_TRACE

LIBDIRS += -DNUOPC

# debug for gfortran
# -- INCDIRS += -g -O0 -fbacktrace 


# debug for intel
# -- INCDIRS += -g -traceback 

LIBSS   :=  -lnetcdf  -lnetcdff -lhdf5 -lhdf5_fortran

FLIBS   := $(LIBDIRS)  $(LIBSS)


########################################################################
# Compiler flags for Linux operating system on 64bit x86 CPU
#
ifeq ($(MACHINE)-$(OS),x86_64-linux-gnu)
#
# ***NOTE*** User must select between various Linux setups
#            by commenting/uncommenting the appropriate compiler
#
#compiler=gnu
#compiler=g95
#compiler=intel
#compiler=intel-ND
#compiler=intel-lonestar
#compiler=cray_xt3
#compiler=cray_xt4
#compiler=cray_xt5
#compiler=pgi
#compiler=pgi-ranger
#compiler=diamond
#compiler=kraken
#compiler=utils
#
#
# Compiler Flags for gfortran and gcc
ifeq ($(compiler),gnu)
  PPFC		:=  gfortran
  FC		:=  gfortran
  PFC		:=  mpif90
  FFLAGS1	:=  $(INCDIRS) -O2 -mcmodel=medium -ffixed-line-length-none -march=k8 -m64
  FFLAGS2	:=  $(FFLAGS1)
  FFLAGS3	:=  $(FFLAGS1)
  DA		:=  -DREAL8 -DLINUX -DCSCA
  DP		:=  -DREAL8 -DLINUX -DCSCA -DCMPI -DHAVE_MPI_MOD
  DPRE		:=  -DREAL8 -DLINUX
  IMODS 	:=  -I
  CC		:= gcc
  CCBE		:= $(CC)
  CFLAGS	:= $(INCDIRS) -O2 -mcmodel=medium -DLINUX -march=k8 -m64
  CLIBS	:=
  LIBS		:=
  MSGLIBS	:=
  ifeq ($(NETCDF),enable)
     ifeq ($(MACHINENAME),blueridge)
        # FLIBS       := $(FLIBS) -L$(HDF5HOME) -lhdf5  
        NETCDFHOME    :=/usr
        FFLAGS1       :=$(FFLAGS1) -I/usr/lib64/gfortran/modules
        FFLAGS2       :=$(FFLAGS1)
        FFLAGS3       :=$(FFLAGS1)
        # NETCDFHOME  :=/shared/apps/RHEL-5/x86_64/NetCDF/netcdf-4.1.1-gcc4.1-ifort
        # NETCDFHOME  :=/shared/apps/RHEL-5/x86_64/NetCDF/netcdf-4.1.2-gcc4.1-ifort
        FLIBS          :=$(FLIBS) -L/usr/lib64 -lnetcdff
     else
        FLIBS          := $(FLIBS) -L$(HDF5HOME) -lhdf5 -lhdf5_fortran
     endif
  endif
  $(warning (INFO) Corresponding compilers and flags found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
endif
#
ifeq ($(compiler),g95)
  PPFC		:=  g95
  FC		:=  g95
  PFC		:=  mpif90
  FFLAGS1	:=  $(INCDIRS) -O3 -mcmodel=medium -fstatic -ffixed-line-length-132
  FFLAGS2	:=  $(FFLAGS1)
  FFLAGS3	:=  $(FFLAGS1)
  DA		:=  -DREAL8 -DLINUX -DCSCA
  DP		:=  -DREAL8 -DLINUX -DCSCA -DCMPI
  DPRE		:=  -DREAL8 -DLINUX
  IMODS 	:=  -I
  CC		:= gcc
  CCBE		:= $(CC)
  CFLAGS	:= $(INCDIRS) -O2 -mcmodel=medium -DLINUX
  CLIBS	:=
  FLIBS		:=
  MSGLIBS	:=
  $(warning (INFO) Corresponding compilers and flags found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
endif
#
# jgf45.12 These flags work on the UNC Topsail Cluster.
# jgf: The -i-dynamic flag defers the inclusion of the library with 
# feupdateenv until run time, thus avoiding the error message:
# "feupdateenv is not implemented and will always fail"
ifeq ($(compiler),intel)
  PPFC          :=  ifort
  FC            :=  ifort
  PFC           :=  mpiifort
#  PFC           :=  mpif90
  FFLAGS1       :=  $(INCDIRS) -w -O3 -FI -assume byterecl -132 -i-dynamic -assume buffered_io
  FFLAGS2       :=  $(FFLAGS1)
  FFLAGS3       :=  $(FFLAGS1)
  DA            :=  -DREAL8 -DLINUX -DCSCA
  DP            :=  -DREAL8 -DLINUX -DCSCA -DCMPI
  DPRE          :=  -DREAL8 -DLINUX
  ifeq ($(SWAN),enable)
     DPRE          := $(DPRE) -DADCSWAN
  endif
  IMODS         :=  -I
  CC            := icc
  CCBE          := $(CC)
  CFLAGS        := $(INCDIRS) -O3 -m64 -mcmodel=medium -DLINUX
  CLIBS         :=
  FLIBS         :=
  MSGLIBS       :=
  NETCDFHOME     :=/apps/netcdf/4.3.0-intel
  HDF5HOME       :=/apps/hdf5/1.8.14-intel/lib/
  ifeq ($(NETCDF),enable)
        FLIBS          := $(FLIBS) -lnetcdff -L$(HDF5HOME) -lhdf5_fortran -lhdf5_hl -lhdf5
  endif
  $(warning (INFO) Corresponding machine found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
endif
#
# Corbitt 120322:  These flags work on the Notre Dame Athos & Zas
ifeq ($(compiler),intel-ND)
  PPFC            :=  ifort
  FC            :=  ifort
  PFC           :=  mpif90
  FFLAGS1       :=  $(INCDIRS) -w -O3 -assume byterecl -132 -i-dynamic -assume buffered_io
  ifeq ($(DEBUG),full)
     FFLAGS1       :=  $(INCDIRS) -g -O0 -traceback -debug -check all -i-dynamic -FI -assume byterecl -132 -DALL_TRACE -DFULL_STACK -DFLUSH_MESSAGES
  endif
  FFLAGS2       :=  $(FFLAGS1)
  FFLAGS3       :=  $(FFLAGS1)
  DA            :=  -DREAL8 -DLINUX -DCSCA
  DP            :=  -DREAL8 -DLINUX -DCSCA -DCMPI -DPOWELL
  DPRE          :=  -DREAL8 -DLINUX -DADCSWAN
  ifeq ($(SWAN),enable)
     DPRE          := $(DPRE) -DADCSWAN
  endif
  IMODS         :=  -I
  CC            := icc
  CCBE          := $(CC)
  CFLAGS        := $(INCDIRS) -O3 -m64 -mcmodel=medium -DLINUX
  FLIBS          :=
  ifeq ($(DEBUG),full)
     CFLAGS        := $(INCDIRS) -g -O0 -march=k8 -m64 -mcmodel=medium -DLINUX
  endif
  ifeq ($(NETCDF),enable)
     HDF5HOME=/afs/crc.nd.edu/x86_64_linux/hdf/hdf5-1.8.6-linux-x86_64-static/lib
     FLIBS      := $(FLIBS) -lnetcdff -L$(HDF5HOME) 
  endif   
  CLIBS         :=
  MSGLIBS       :=
  $(warning (INFO) Corresponding machine found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
  NETCDFHOME=/afs/crc.nd.edu/x86_64_linux/netcdf/rhel6/4.1.3/intel-12.0/
  #NETCDFHOME=/afs/crc.nd.edu/x86_64_linux/scilib/netcdf/4.1.2/intel-12.0/inst
endif
#
# sb46.50.02 These flags work on the UT Austin Lonstar cluster.
ifeq ($(compiler),intel-lonestar)
  PPFC            :=  ifort
  FC            :=  ifort
  PFC           :=  mpif90
  FFLAGS1       :=  $(INCDIRS) -O3 -xT -132 -i-dynamic
  ifeq ($(DEBUG),full)
     FFLAGS1       :=  $(INCDIRS) -g -O0 -traceback -debug all -ftrapuv -fpe0 -check all -i-dynamic -FI -assume byterecl -132 -DALL_TRACE -DFULL_STACK -DFLUSH_MESSAGES
  endif
  FFLAGS2       :=  $(FFLAGS1)
  FFLAGS3       :=  $(FFLAGS1)
  DA            :=  -DREAL8 -DLINUX -DCSCA
  DP            :=  -DREAL8 -DLINUX -DCSCA -DCMPI
  DPRE          :=  -DREAL8 -DLINUX
  ifeq ($(SWAN),enable)
     DPRE          := $(DPRE) -DADCSWAN
  endif
  IMODS         :=  -I
  CC            := icc
  CCBE		:= $(CC)
  CFLAGS        := $(INCDIRS) -O3 -xT
  ifeq ($(DEBUG),full)
     CFLAGS        := $(INCDIRS) -g -O0 
  endif
  CLIBS         :=
  FLIBS          :=
  MSGLIBS       :=
  $(warning (INFO) Corresponding machine found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
  NETCDFHOME=/opt/apps/intel11_1/netcdf/4.2.1.1/
  HDF5HOME=/opt/apps/intel11_1/hdf5/1.8.8/
  TACC_NETCDF_INC=${NETCDFHOME}include
  TACC_NETCDF_LIB=${NETCDFHOME}lib
  TACC_HDF5_LIB=${HDF5HOME}lib
  ifeq ($(NETCDF),enable)
     FLIBS      := $(FLIBS) -I${TACC_NETCDF_INC} -L${TACC_NETCDF_LIB} -lnetcdf -lnetcdff -L${TACC_HDF5_LIB} -lhdf5 -lhdf5 -lz -lm
  endif
endif
#
# Cray-XT3 using standard compilers, from vjp; added by jgf46.00
ifeq ($(compiler),cray_xt3)
  PPFC	        :=  pgf90
  FC	        :=  ftn
  PFC	        :=  ftn
  CC		:=  pgcc
  CCBE		:=  cc
  FFLAGS1	:=  $(INCDIRS) -fastsse -Mextend
  FFLAGS2	:=  $(FFLAGS1)
  FFLAGS3	:=  $(FFLAGS1) -r8 -Mr8 -Mr8intrinsics
  DA  	        :=  -DREAL8 -DLINUX -DCSCA
  DP  	        :=  -DREAL8 -DLINUX -DCMPI -DHAVE_MPI_MOD -DCSCA -DDEBUG_WARN_ELEV
#  DP  	        :=  -DREAL8 -DLINUX -DCMPI -DHAVE_MPI_MOD -DCSCA
  DPRE	        :=  -DREAL8 -DLINUX
  CFLAGS	:=  -c89 $(INCDIRS) -DLINUX
  IMODS		:=  -module
  FLIBS  	:=
  MSGLIBS	:=
# When compiling with netCDF support, the HDF5 libraries must also
# be linked in, so the user must specify HDF5HOME on the command line.
# jgf20101102: on Sapphire,
#              NETCDFHOME=/usr/local/usp/PETtools/CE/pkgs/netcdf-4.1.1-serial
#              HDF5HOME=${PET_HOME}/pkgs/hdf5-1.8.5-serial/lib
  ifeq ($(NETCDF),enable)
     FLIBS          := $(FLIBS) -L$(HDF5HOME) -lhdf5_fortran -lhdf5_hl -lhdf5 -lz
  endif
  BACKEND_EXEC  := metis_be adcprep_be
  $(warning (INFO) Corresponding machine found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
endif
#
# Cray-XT4 (e.g. Jade@ERDC) using standard compilers, added by jgf48.4607
ifeq ($(compiler),cray_xt4)
  PPFC	        :=  pgf90
  FC	        :=  ftn
  PFC	        :=  ftn
  CC		:=  pgcc
  CCBE		:=  cc
  FFLAGS1	:=  $(INCDIRS) -Mextend -Minform,inform -O2 -fastsse 
  ifeq ($(DEBUG),full)
     FFLAGS1	:=  $(INCDIRS) -Mextend -g -O0 -traceback -Mbounds -Mchkfpstk -Mchkptr -Mchkstk -DALL_TRACE -DFLUSH_MESSAGES -DFULL_STACK
  endif
  FFLAGS2	:=  $(FFLAGS1)
  FFLAGS3	:=  $(FFLAGS1) -r8 -Mr8 -Mr8intrinsics
  DA  	        :=  -DREAL8 -DLINUX -DCSCA
  DP  	        :=  -DREAL8 -DLINUX -DCMPI -DHAVE_MPI_MOD -DCSCA
  DPRE	        :=  -DREAL8 -DLINUX
  ifeq ($(SWAN),enable)
     DPRE	        :=  -DREAL8 -DLINUX -DADCSWAN
  endif
  CFLAGS	:=  $(INCDIRS) -DLINUX
  ifeq ($(DEBUG),full)
     CFLAGS	:=  $(INCDIRS) -DLINUX -g -O0
  endif
  IMODS		:=  -module
  FLIBS  	:=
# When compiling with netCDF support, the HDF5 libraries must also
# be linked in, so the user must specify HDF5HOME on the command line.
# On Jade, HDF5 was compiled with szip compression, so this library is
# required as well.
# jgf20101102: on Jade, NETCDFHOME=/usr/local/usp/PETtools/CE/pkgs/netcdf-4.0.1-serial
# jgf20101102: on Jade, HDF5HOME=${PET_HOME}/pkgs/hdf5-1.8.4-serial/lib
# jgf20101103: on Jade, SZIPHOME=/usr/local/usp/PETtools/CE/pkgs/szip-2.1/lib
# jgf20110728: on Garnet, NETCDFHOME=/opt/cray/netcdf/4.1.1.0/netcdf-pgi
  ifeq ($(NETCDF),enable)
     FLIBS          := $(FLIBS) -L$(HDF5HOME) -L$(SZIPHOME) -lhdf5_fortran -lhdf5_hl -lhdf5 -lsz -lz
  endif
  MSGLIBS	:=
  BACKEND_EXEC  := metis_be adcprep_be
  $(warning (INFO) Corresponding machine found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
endif
#
# Cray-XT5 (e.g. einstein@NAVO) using standard compilers, added by jgf49.07
ifeq ($(compiler),cray_xt5)
  PPFC	        :=  ftn
  FC	        :=  ftn
  PFC	        :=  ftn
  CC		:=  cc
  CCBE		:=  cc
  FFLAGS1	:=  $(INCDIRS) -Mextend -Minform,inform -O2 -fastsse
#  FFLAGS1	:=  $(INCDIRS) -Mextend -g -O0 -traceback
  FFLAGS2	:=  $(FFLAGS1)
  FFLAGS3	:=  $(FFLAGS1) -r8 -Mr8 -Mr8intrinsics
  DA  	        :=  -DREAL8 -DLINUX -DCSCA
  DP  	        :=  -DREAL8 -DLINUX -DCMPI -DHAVE_MPI_MOD -DCSCA
  DPRE	        :=  -DREAL8 -DLINUX
  CFLAGS	:=  -c89 $(INCDIRS) -DLINUX
  IMODS		:=  -module
  FLIBS  	:=
# When compiling with netCDF support, the HDF5 libraries must also
# be linked in, so the user must specify HDF5HOME on the command line.
# jgf20090518: on Jade, NETCDFHOME=/usr/local/usp/PETtools/CE/pkgs/netcdf-4.0
# jgf20090518: on Jade, HDF5HOME=${PET_HOME}/pkgs/hdf5-1.8.2/lib
  ifeq ($(NETCDF),enable)
     FLIBS          := $(FLIBS) -L$(HDF5HOME) -lhdf5 -lhdf5_fortran
  endif
  MSGLIBS	:=
  BACKEND_EXEC  := metis_be adcprep_be
  $(warning (INFO) Corresponding machine found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
endif
#

# Portland Group
ifeq ($(compiler),pgi)
  PPFC		:=  pgf90
  FC		:=  pgf90
  PFC		:=  mpif90
  FFLAGS1	:=  $(INCDIRS) -fastsse -mcmodel=medium -Mextend
  FFLAGS2	:=  $(FFLAGS1)
  FFLAGS3	:=  $(FFLAGS1)
  DA		:=  -DREAL8 -DLINUX -DCSCA
  DP		:=  -DREAL8 -DLINUX -DCSCA -DCMPI
  DPRE		:=  -DREAL8 -DLINUX
  IMODS		:=  -I
  CC		:= gcc
  CCBE          := $(CC)
  CFLAGS	:= $(INCDIRS) -O2 -mcmodel=medium -DLINUX
  CLIBS		:=
  FLIBS  	:=
  MSGLIBS	:=
  $(warning (INFO) Corresponding machine found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
endif
#
# Utility Server (e.g. US@ERDC) using standard compilers, added by jgf48.4607
ifeq ($(compiler),utils)
  PPFC          :=  pgf90
  FC            :=  pgf90
  PFC           :=  mpif90
  CC            :=  pgcc
  CCBE          :=  pgcc
  FFLAGS1       :=  $(INCDIRS) -Mextend -Minform,inform -O2 -fastsse
  ifeq ($(DEBUG),full)
     FFLAGS1    :=  $(INCDIRS) -Mextend -g -O0 -traceback -Mbounds -Mchkfpstk -Mchkptr -Mchkstk -DALL_TRACE -DFLUSH_MESSAGES -DFULL_STACK
  endif
  FFLAGS2       :=  $(FFLAGS1)
  FFLAGS3       :=  $(FFLAGS1) -r8 -Mr8 -Mr8intrinsics
  DA            :=  -DREAL8 -DLINUX -DCSCA
  DP            :=  -DREAL8 -DLINUX -DCMPI -DHAVE_MPI_MOD -DCSCA
  DPRE          :=  -DREAL8 -DLINUX
  ifeq ($(SWAN),enable)
     DPRE               :=  -DREAL8 -DLINUX -DADCSWAN
  endif
  CFLAGS        :=  -c89 $(INCDIRS) -DLINUX
  ifeq ($(DEBUG),full)
     CFLAGS     :=  -c89 $(INCDIRS) -DLINUX -g -O0
  endif
  IMODS         :=  -module
  FLIBS         :=
# When compiling with netCDF support, the HDF5 libraries must also
# be linked in, so the user must specify HDF5HOME on the command line.
# On Jade, HDF5 was compiled with szip compression, so this library is
# required as well.
# jgf20101102: on Jade, NETCDFHOME=/usr/local/usp/PETtools/CE/pkgs/netcdf-4.0.1-serial
# jgf20101102: on Jade, HDF5HOME=${PET_HOME}/pkgs/hdf5-1.8.4-serial/lib
# jgf20101103: on Jade, SZIPHOME=/usr/local/usp/PETtools/CE/pkgs/szip-2.1/lib
# jgf20110728: on Garnet, NETCDFHOME=/opt/cray/netcdf/4.1.1.0/netcdf-pgi
  ifeq ($(NETCDF),enable)
     FLIBS          := $(FLIBS) -L$(HDF5HOME) -L$(SZIPHOME) -lhdf5_fortran -lhdf5_hl -lhdf5 -lsz -lz
  endif
  MSGLIBS       :=
  BACKEND_EXEC  := metis_be adcprep_be
  $(warning (INFO) Corresponding machine found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
endif
# Portland Group on TU Ranger (AMD Opteron 8356, Barcelona Core)  Seizo
ifeq ($(compiler),pgi-ranger)
  PPFC          :=  pgf95
  FC            :=  pgf95
  PFC           :=  mpif90
  FFLAGS1       :=  $(INCDIRS) -fast -tp barcelona-64 -Mextend
  ifeq ($(DEBUG),full)
     FFLAGS1	:=  $(INCDIRS) -Minform,inform -Mextend -g -O0 -traceback -DNETCDF_DEBUG -Mbounds -Mchkfpstk -Mchkptr -Mchkstk -DALL_TRACE -DFLUSH_MESSAGES -DFULL_STACK
  endif
  FFLAGS2       :=  $(FFLAGS1)
  FFLAGS3       :=  $(FFLAGS1)
  DA            :=  -DREAL8 -DLINUX -DCSCA
  DP            :=  -DREAL8 -DLINUX -DCSCA -DCMPI
  DPRE          :=  -DREAL8 -DLINUX
  IMODS         :=  -I
  CC            := gcc
  CCBE          := $(CC)
  CFLAGS        := $(INCDIRS) -DLINUX
  CLIBS         :=
  LIBS          :=
  MSGLIBS       :=
  $(warning (INFO) Corresponding machine found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
endif
#
# ERDC Diamond
ifeq ($(compiler),diamond)
  PPFC          :=  ifort
  FC            :=  ifort
  PFC           :=  ifort
#  FFLAGS1       :=  $(INCDIRS) -O3 -xT -132
  FFLAGS1       := -O3 -132 -xSSSE3
  ifeq ($(DEBUG),full)
     FFLAGS1	:=  $(INCDIRS) -g -O0 -debug -fpe0 -132 -traceback -check all -DALL_TRACE -DFLUSH_MESSAGES -DFULL_STACK
  endif
  FFLAGS2       :=  $(FFLAGS1)
  FFLAGS3       :=  $(FFLAGS1)
  DA            :=  -DREAL8 -DLINUX -DCSCA
  DP            :=  -DREAL8 -DLINUX -DCSCA -DCMPI -DHAVE_MPI_MOD
  DPRE          :=  -DREAL8 -DLINUX
  ifeq ($(SWAN),enable)
     DPRE          :=  -DREAL8 -DLINUX -DADCSWAN
  endif
  IMODS         :=  -I
  CC            := icc
  CCBE          := $(CC)
#  CFLAGS        := $(INCDIRS) -O3 -xT
  CFLAGS        := $(INCDIRS) -O3 -xSSSE3
  ifeq ($(DEBUG),full)
     CFLAGS        := $(INCDIRS) -g -O0
  endif
  CLIBS         :=
  LIBS          :=
# When compiling with netCDF support, the HDF5 libraries must also
# be linked in, so the user must specify HDF5HOME on the command line.
# jgf20101103: on Diamond, NETCDFHOME=/usr/local/usp/PETtools/CE/pkgs/netcdf-4.0.1-serial
# jgf20101103: on Diamond, HDF5HOME=${PET_HOME}/pkgs/hdf5-1.8.4-serial/lib
  ifeq ($(NETCDF),enable)
     NETCDFHOME     :=/usr/local/usp/PETtools/CE/pkgs/netcdf-4.2.1.1-intel-serial
     FFLAGS1        :=$(FFLAGS1) -I/usr/local/usp/PETtools/CE/pkgs/netcdf-4.2.1.1-intel-serial/include
     FFLAGS2        :=$(FFLAGS1)
     FFLAGS3        :=$(FFLAGS1)
     #HDF5HOME       :=/usr/local/usp/PETtools/CE/pkgs/hdf5-1.8.8-serial/lib
     FLIBS          := $(FLIBS) -L/usr/local/usp/PETtools/CE/pkgs/netcdf-4.2.1.1-intel-serial/lib -lnetcdff -lnetcdf   #-L$(HDF5HOME) -lhdf5_hl -lhdf5 -lhdf5_fortran -lz
  endif
  MSGLIBS       := -lmpi
  $(warning (INFO) Corresponding machine found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
endif
#
# Cray-XE6 (e.g., Garnet at ERDC) using standard compilers, added by jgf50.29
ifeq ($(compiler),garnet)
  PPFC	        :=  ftn
  FC	        :=  ftn
  PFC	        :=  ftn
  CC		:=  pgcc
  CCBE		:=  cc
  FFLAGS1	:=  $(INCDIRS) -Mextend -Minform,inform -tp=shanghai-64 -fast
  ifeq ($(DEBUG),full)
     FFLAGS1	:=  $(INCDIRS) -Mextend -g -O0 -traceback -Mbounds -Mchkfpstk -Mchkptr -Mchkstk -DALL_TRACE -DFLUSH_MESSAGES -DFULL_STACK -tp=shanghai-64
  endif
  FFLAGS2	:=  $(FFLAGS1)
  FFLAGS3	:=  $(FFLAGS1) #-r8 -Mr8 -Mr8intrinsics
  DA  	        :=  -DREAL8 -DLINUX -DCSCA
  DP  	        :=  -DREAL8 -DLINUX -DCMPI -DHAVE_MPI_MOD -DCSCA
  DPRE	        :=  -DREAL8 -DLINUX
  ifeq ($(SWAN),enable)
     DPRE	        :=  -DREAL8 -DLINUX -DADCSWAN
  endif
  CFLAGS	:=  $(INCDIRS) -DLINUX
  ifeq ($(DEBUG),full)
     CFLAGS	:=  $(INCDIRS) -DLINUX -g -O0
  endif
  IMODS		:=  -module
  FLIBS         := 
# jgf20110728: on Garnet, NETCDFHOME=/opt/cray/netcdf/4.1.1.0/netcdf-pgi
# jgf20110815: on Garnet, HDF5HOME=/opt/cray/hdf5/default/hdf5-pgi
# jgf20130815: on Garnet, load module cray-netcdf, with the path to the
#              installation being /opt/cray/netcdf/4.3.0
  ifeq ($(NETCDF),enable)
     FLIBS          := $(FLIBS) -lnetcdff
  endif
  FLIBS  	:= $(FLIBS)
  MSGLIBS	:=
  BACKEND_EXEC  := metis_be adcprep_be
  $(warning (INFO) Corresponding machine found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
endif
#
#
# Tennessee's Kraken
ifeq ($(compiler),kraken)
  PPFC          :=  ftn
  FC            :=  ftn
  PFC           :=  ftn
  FFLAGS1       :=  $(INCDIRS) -O3 -static  -132
  FFLAGS2       :=  $(FFLAGS1)
  FFLAGS3       :=  $(FFLAGS1)
  DA            :=  -DREAL8 -DLINUX -DCSCA -DPOWELL
  DP            :=  -DREAL8 -DLINUX -DCSCA -DCMPI -DPOWELL
  DPRE          :=  -DREAL8 -DLINUX -DADCSWAN
  IMODS         :=  -I
  CC            := cc
  CCBE          := $(CC)
  CFLAGS        := $(INCDIRS)  -DLINUX
  CLIBS         :=
  LIBS          :=
  MSGLIBS       :=
  $(warning (INFO) Corresponding machine found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
endif
endif
#$(MACHINE)
########################################################################
# Compiler flags for Linux operating system on 32bit x86 CPU
#
ifeq ($(MACHINE)-$(OS),i686-linux-gnu)
#
# ***NOTE*** User must select between various Linux setups
#            by commenting/uncommenting the appropriate compiler or
#            by specifying on the make command line, e.g.,
#            make all compiler=gnu
#
#compiler=gnu
#compiler=intel
#compiler=pgi
#compiler=gfortran
#
# Portland Group
ifeq ($(compiler),pgi)
  PPFC          :=  pgf90
  FC	        :=  pgf90
  PFC	        :=  mpif90
  FFLAGS1	:=  $(INCDIRS) -O2
  FFLAGS2	:=  $(FFLAGS1)
  FFLAGS3	:=  $(FFLAGS1)
  DA  	        :=  -DREAL8 -DLINUX -DCSCA
  DP  	        :=  -DREAL8 -DLINUX -DCSCA -DCMPI
  DPRE	        :=  -DREAL8 -DLINUX
  IMODS 	:=  -I
  CC            := pgcc
  CCBE          := $(CC)
  CFLAGS       := $(INCDIRS) -O2 -DLINUX
  CLIBS         :=
  FLIBS  	:=
  MSGLIBS	:=
  $(warning (INFO) Corresponding machine found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
endif
#
# Intel
ifeq ($(compiler),intel)
  PPFC	        :=  ifort -w
  FC	        :=  ifort -w
  PFC	        :=  mpif90
  OPTLVL        := -O2
  ifeq ($(ADC_DEBUG),yes)
    OPTLVL        := -g
  endif
  FFLAGS1	:=  $(INCDIRS) $(OPTLVL) -extend_source -Vaxlib -assume byterecl
  FFLAGS2	:=  $(FFLAGS1)
  FFLAGS3	:=  $(FFLAGS1)
  DA  	        :=  -DREAL8 -DLINUX -DCSCA
  DP  	        :=  -DREAL8 -DLINUX -DCSCA -DCMPI -DHAVE_MPI_MOD
  DPRE	        :=  -DREAL8 -DLINUX
  IMODS 	:=  -I
  CC            := icc
  CCBE          := $(CC)
  CFLAGS        := $(INCDIRS) $(OPTLVL) -DLINUX
  CLIBS         :=
  FLIBS  	:=
  MSGLIBS	:=
  $(warning (INFO) Corresponding machine found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
endif
#
# g95
ifeq ($(compiler),gnu)
  PPFC		:=  g95
  FC		:=  g95
  PFC		:=  mpif90
  FFLAGS1	:=  $(INCDIRS) -O2 -ffixed-line-length-132
  ifeq ($(DEBUG),full)
     FFLAGS1	:=  $(INCDIRS) -g -O0 -ffixed-line-length-132 -ftrace=full -fbounds-check -DALL_TRACE -DFLUSH_MESSAGES -DFULL_STACK -DWRITER_DEBUG -DDEBUG_HOLLAND
      # g95 environment variables to set for enhanced debugging:
      # G95_UNBUFFERED_ALL, G95_ABORT, G95_FPU_DENORMAL, G95_FPU_INVALID,
      # G95_FPU_ZERODIV, G95_FPU_OVERFLOW, G95_FPU_UNDERFLOW,
      # G95_FPU_EXCEPTIONS
  endif
  ifeq ($(DEBUG),netcdf)
     FFLAGS1	:=  $(INCDIRS) -g -O0 -ffixed-line-length-132 -ftrace=full -fbounds-check -DNETCDF_TRACE -DFLUSH_MESSAGES -DFULL_STACK
  endif
  ifeq ($(DEBUG),valgrind)
     FFLAGS1	:=  $(INCDIRS) -g -O0 -ffixed-line-length-132 
  endif
  ifeq ($(SWAN),enable)
     FFLAGS1    :=  $(FFLAGS1) -freal-loops
  endif
  FFLAGS2	:=  $(FFLAGS1)
  FFLAGS3	:=  $(FFLAGS1)
  DA		:=  -DREAL8 -DLINUX -DCSCA
  DP		:=  -DREAL8 -DLINUX -DCSCA -DCMPI
  DPRE		:=  -DREAL8 -DLINUX
  ifeq ($(SWAN),enable)
     DPRE               :=  -DREAL8 -DLINUX -DADCSWAN
  endif
  IMODS 	:=  -I
  CC		:= gcc
  CCBE		:= $(CC)
  CFLAGS	:= $(INCDIRS) -O2 -DLINUX
  ifeq ($(DEBUG),full)
     CFLAGS     := $(INCDIRS) -g -O0 -DLINUX
  endif
  CLIBS	:=
  FLIBS		:=
  ifeq ($(NETCDF),enable)
     ifeq ($(MACHINENAME),jason-desktop)
        NETCDFHOME := /usr/local
        HDF5HOME   := /usr/local/hdf5/hdf5-1.8.8/hdf5
     endif
     FLIBS          := $(FLIBS) -L$(HDF5HOME)/lib -L$(NETCDFHOME) -lnetcdf -lhdf5_hl -lhdf5 -lhdf5_fortran -lz
  endif
  MSGLIBS	:=
  $(warning (INFO) Corresponding machine found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
endif
#
# gfortran
ifeq ($(compiler),gfortran)
  ifeq ($(MACHINENAME),jason-desktop)
     XDMFPATH    := /home/jason/projects/XDMF/Code/latestCode
     XDMFLIBPATH := /home/jason/projects/XDMF/Code/testLatest
  endif
  PPFC		:=  gfortran
  FC		:=  gfortran
  PFC		:=  mpif90
  FFLAGS1	:=  $(INCDIRS) -O2 -ffixed-line-length-none 
  ifeq ($(PROFILE),enable)
    FFLAGS1	:=  $(INCDIRS) -pg -O0 -fprofile-arcs -ftest-coverage -ffixed-line-length-none 
  endif
  ifeq ($(DEBUG),full)
    FFLAGS1	:=  $(INCDIRS) -g -O0 -ffixed-line-length-none -fbacktrace -fbounds-check -ffpe-trap=zero,invalid,overflow,denormal -DALL_TRACE -DFLUSH_MESSAGES -DFULL_STACK -DDEBUG_HOLLAND
  endif
  ifeq ($(DEBUG),full-not-fpe)
    FFLAGS1	:=  $(INCDIRS) -g -O0 -ffixed-line-length-none -fbacktrace -fbounds-check -DALL_TRACE -DFLUSH_MESSAGES -DFULL_STACK -DDEBUG_HOLLAND
  endif
  ifneq ($(MACHINENAME),jason-desktop)
     FFLAGS1 := $(FFLAGS1) -fno-underscoring
  endif
  FFLAGS2	:=  $(FFLAGS1)
  FFLAGS3	:=  $(FFLAGS1)
  DA		:=  -DREAL8 -DLINUX -DCSCA
  DP		:=  -DREAL8 -DLINUX -DCSCA -DCMPI
  DPRE		:=  -DREAL8 -DLINUX
  ifeq ($(SWAN),enable)
     DPRE               :=  -DREAL8 -DLINUX -DADCSWAN
  endif
  FLIBS         := 
  ifeq ($(NETCDF),enable)
     ifeq ($(MACHINENAME),jason-desktop)
        NETCDFHOME := /usr
     endif
     FLIBS      := $(FLIBS) -lnetcdff 
  endif
  IMODS 	:=  -I
  CC		:= gcc
  CCBE		:= $(CC)
  CFLAGS	:= $(INCDIRS) -O2 -DLINUX
  ifeq ($(DEBUG),full)
     CFLAGS     := $(INCDIRS) -g -O0 -DLINUX
  endif
  CLIBS	:=
  MSGLIBS	:=
  $(warning (INFO) Corresponding machine found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
endif

endif

########################################################################
# Compiler flags for Linux operating system on 64bit Intel Itanium CPU
#
ifneq (,$(findstring ia64-linux,$(MACHINE)-$(OS)))
#
# ***NOTE*** User must select between various Linux setups
#            by commenting/uncommenting the appropriate architecture
#
arch=pc
#arch=altix
#
#
ifeq ($(arch),pc)
  PPFC          := /opt/intel/compiler60/ia64/bin/efc
  FC            := /opt/intel/compiler60/ia64/bin/efc
  PFC           :=  mpif90
  FFLAGS1       := $(INCDIRS) -132 -IPF_fp_speculationfast -r8
  FFLAGS2       := $(FFLAGS1)
  FFLAGS3       := $(FFLAGS1)
  DA            :=  -DREAL8 -DCSCA
  DP            :=  -DREAL8 -DCSCA -DCMPI
  DPRE          :=  -DREAL8 -DLINUX
  IMODS         :=  -I
  CC            :=  gcc
  CCBE          :=  $(CC)
  CFLAGS        :=  $(INCDIRS) -DLINUX -O2
  CLIBS         :=
  FLIBS          :=
  MSGLIBS       :=
  $(warning (INFO) Corresponding machine found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
endif
#
# SGI Altix; Created by Brett Estrade, added by jgf45.12
ifeq ($(arch),altix)
  PPFC            := ifort
  FC              := ifort
  PFC             := ifort
  FFLAGS1	  := $(INCDIRS) -O3 -tpp2
  FFLAGS2	  := $(FFLAGS1)
  FFLAGS3	  := $(FFLAGS1)
  DA	          :=  -DREAL8 -DCSCA
  DP	          :=  -DREAL8 -DCSCA -DCMPI
  DPRE	          :=  -DREAL8 -DLINUX
  IMODS   	  :=  -I
  CC              :=  gcc
  CCBE            :=  $(CC)
  CFLAGS          :=  $(INCDIRS) -DLINUX -O2
  FLIBS		  :=
  MSGLIBS	  := -lmpi
  $(warning (INFO) Corresponding machine found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
endif
endif
########################################################################
# IBM SP - AIX operating system on IBM RS/6000 hardware
# Contrib. TJC Jan 2006
ifneq (,$(findstring rs6000-aix,$(MACHINE)-$(OS)))
   export OBJECT_MODE=64
   PPFC          := xlf90_r
   FC            := xlf90_r
   PFC           := mpxlf90_r
   FFLAGS1       := $(INCDIRS) -q64 -w -O2 -qfixed=132 -qarch=auto -qcache=auto
   FFLAGS2       := $(FFLAGS1)
   FFLAGS3       := $(FFLAGS1)
   DA            := -WF,"-DREAL8,-DIBM,-DCSCA"
   DP            := -tF -WF,"-DREAL8,-DIBM,-DCSCA,-DCMPI"
   DPRE          := -tF -WF,"-DREAL8,-DIBM"
   IMODS         := -I
   CC            := mpcc_r
   CCBE          :=  $(CC)
   CFLAGS        := -q64 -I. -O2 -DIBM
   LDFLAGS       := -q64
   FLIBS          :=
   MSGLIBS       := -lm
  $(warning (INFO) Corresponding machine found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
endif
########################################################################

########################################################################
# IBM SP - AIX operating system on 32bit PowerPC CPU
# From vjp; added by jgf46.00
# gcc by bde for IBM p5: bluedawg.loni.org; added by jgf46.04
#
# jgf20110303: Added capability to specify ibm=p5 or ibm=p6 on make
# command line to get the appropriate compiler flags. It also possible
# to uncomment either ibm=p5 or ibm=p6 below and not specify it on the
# make command line. The ibm p6 compiler flags were inserted based on
# feedback from Yuji Funakoshi at NOAA CSDL.
#
ifneq (,$(findstring powerpc-aix,$(MACHINE)-$(OS)))
#IBM=p5
#IBM=p6
ifeq ($(IBM),p5)
  PPFC          := xlf90_r -q64
  FC            := xlf90_r -q64
  PFC           := mpxlf90_r -q64
  CC            := mpcc_r  -q64
  CCBE          :=  $(CC)
  FFLAGS0       := $(INCDIRS) -w -qfixed=132 -qarch=auto -qcache=auto
  FFLAGS1       := $(FFLAGS0) -O2
  FFLAGS2       := $(FFLAGS0) -qhot -qstrict
  FFLAGS3       := $(FFLAGS0) -O3 -qinitauto
  DA            := -WF,"-DREAL8,-DIBM,-DCSCA"
  DP            := -tF -WF,"-DREAL8,-DIBM,-DCSCA,-DCMPI"
  DPRE          := -tF -WF,"-DREAL8,-DIBM"
  IMODS         := -I
  CFLAGS        := $(INCDIRS) -O2 -DIBM
  ARFLAGS	:= -X64 rv
  FLIBS          :=
  MSGLIBS       := -lm

  $(warning (INFO) Corresponding machine found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
endif
ifeq ($(IBM),p6)
   export OBJECT_MODE=64
   PPFC          := xlf90_r
   FC            := xlf90_r
   PFC           := mpxlf90_r
   FFLAGS1       := $(INCDIRS) -q64 -w -O2 -qfixed=132 -qarch=auto -qcache=auto
   FFLAGS2       := $(FFLAGS1)
   FFLAGS3       := $(FFLAGS1)
   ifeq ($(NETCDF),enable))
      DA            := -WF,"-DREAL8,-DIBM,-DNETCDF,-DCSCA"
      DP            := -tF -WF,"-DREAL8,-DIBM,-DNETCDF,-DCSCA,-DCMPI"
      DPRE          := -tF -WF,"-DREAL8,-DIBM,-DNETCDF"
   else
      DA            := -WF,"-DREAL8,-DIBM,-DCSCA"
      DP            := -tF -WF,"-DREAL8,-DIBM,-DCSCA,-DCMPI"
      DPRE          := -tF -WF,"-DREAL8,-DIBM"
   endif
   IMODS         := -I
   CC            := mpcc_r
   CCBE          :=  $(CC)
   CFLAGS        := -q64 -I. -O2 -DIBM
   LDFLAGS       := -q64
   FLIBS          :=
   MSGLIBS       := -lm
# When compiling with netCDF support, the HDF5 libraries must also
# be linked in, so the user must specify HDF5HOME on the command line.
# yf20110301: on Cirrus/Stratus, NETCDFHOME=/usrx/local/bin/
# yf20110301: on Cirrus/Stratus, HDF5HOME=/usrx/local/hdf5/lib
  HDF5HOME=/usrx/local/hdf5/lib
  ifeq ($(NETCDF),enable)
     FLIBS          := $(FLIBS) -L$(HDF5HOME) -lhdf5_fortran -lhdf5_hl -lhdf5 -lz
  endif
  $(warning (INFO) Corresponding machine found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
endif
endif

########################################################################
# IBM BlueGene/L; Linux operating system on 64bit PowerPC CPU
# Contrib. by Brett Estrade added by jgf45.12
ifeq ($(MACHINE)-$(OS),ppc64-unknown-linux-gnu)
   PPFC          := blrts_xlf90
   FC            := blrts_xlf90
   PFC           := $(FC)
   FFLAGS1       := $(INCDIRS) -qmaxmem=64000 -qfixed -I/bgl/BlueLight/ppcfloor/bglsys/include -O5 -qtune=440 -qcache=440 -qarch=440
   FFLAGS2       := $(FFLAGS1)
   FFLAGS3       := $(FFLAGS1)
   DA            := -WF,"-DREAL8,-DIBM"
   DP            := -tF -WF,"-DREAL8,-DIBM,-DCSCA,-DCMPI"
   DPRE          := -tF -WF,"-DREAL8,-DIBM"
   IMODS         := -I
   CC            := blrts_xlc
   CCBE          :=  $(CC)
   CFLAGS        := -I. -I../Lib -O5
   LDFLAGS       := -q64 -L/bgl/BlueLight/ppcfloor/bglsys/lib -lmsglayer.rts -lrts.rts -ldevices.rts
   FLIBS          :=
   MSGLIBS       := -L/bgl/BlueLight/ppcfloor/bglsys/lib -lmpich.rts -lmsglayer.rts -lrts.rts -ldevices.rts

  $(warning (INFO) Corresponding machine found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
endif
########################################################################



########################################################################
# Sun Solaris
#
#ifneq (,$(findstring sparc-solaris,$(MACHINE)-$(OS)))
#  PPFC	        := f90
#  FC	        := f90
#  PFC	        := mpf90
#  LD	        := $(FC)
#  FFLAGS1	:= -O2
#  FFLAGS2	:= -O2
#  FFLAGS3	:= -O2
#  DA	        := -DREAL8 -DMACHSUN -DCSCA
#  DP	        := -DREAL8 -DMACHSUN -DCSCA -DCMPI
#  DPRE          := -DREAL8 -DMACHSUN
#  IMODS		:= -M
#  CC       	:= cc
#  CCBE          :=  $(CC)
#  CFLAGS   	:= -I. -xO2 -DMACHSUN
#  LIBS     	:=
#  MSGLIBS  	:= -lmpi
#endif

########################################################################
# Sun Solaris Using SUN MPI?

ifneq (,$(findstring sparc-solaris,$(MACHINE)-$(OS)))
  PPFC	        := tmf90
  FC	        := tmf90
  PFC	        := tmf90
  ARCH	        := v8plusa
  LD	        := $(FC)
  FFLAGS1	:= $(INCDIRS) -fast -xO4 -xdepend -fsimple=1 -f -dalign -xtarget=ultra -xarch=$(ARCH)
  FFLAGS2	:= $(FFLAGS1)
  FFLAGS3	:= $(FFLAGS1)
  DA	        := -DREAL8 -DMACHSUN -DCSCA
  DP	        := -DREAL8 -DCMPI -DMACHSUN -DCSCA
  DPRE          := -DREAL8 -DMACHSUN
  IMODS		:= -M
  CC       	:= tmcc
  CCBE          :=  $(CC)
  CFLAGS   	:= $(INCDIRS)
  FLIBS     	:=
  MSGLIBS  	:= -lmpi

  $(warning (INFO) Corresponding machine found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
endif

########################################################################
# Alpha-Linux computers using Compaq (a.k.a. DEC)

ifneq (,$(findstring alphaev6-linux,$(MACHINE)-$(OS)))
  PPFC	        :=  fort
  FC	        :=  fort
  PFC	        :=  mpif90
  FFLAGS1	:=  $(INCDIRS) -fixed -O2
  FFLAGS2	:=  -fixed -O2
  FFLAGS3	:=  -fixed -O2
  DA  	        :=  -DREAL8 -DLINUX -DCSCA
  DP  	        :=  -DREAL8 -DLINUX -DCSCA -DCMPI
  DPRE	        :=  -DREAL8 -DLINUX
  IMODS 	:=  -I
  CC            := ccc
  CCBE          :=  $(CC)
  CFLAGS        := $(INCDIRS) -O2 -DLINUX
  CLIBS         :=
  FLIBS  	:=
  MSGLIBS	:=
  $(warning (INFO) Corresponding machine found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
endif

########################################################################
# Compaq True64 computers

ifneq (,$(findstring alphaev6-osf,$(MACHINE)-$(OS)))
  PPFC		:=  f90
  FC		:=  f90
  PFC		:=  f90
  FFLAGS1	:=  $(INCDIRS) -fixed -O2
  FFLAGS2	:=  $(FFLAGS1)
  FFLAGS3	:=  $(FFLAGS1)
  DA		:=  -DREAL8 -DLINUX -DCSCA
  DP		:=  -DREAL8 -DLINUX -DCSCA -DCMPI
  DPRE		:=  -DREAL8 -DLINUX
  IMODS		:=  -I
  CC		:= cc
  CCBE          :=  $(CC)
  CFLAGS	:= $(INCDIRS) -DLINUX -O2
  CLIBS		:=
  FLIBS		:=
  MSGLIBS	:= -lfmpi -lmpi -lelan
  $(warning (INFO) Corresponding machine found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
endif


########################################################################
# CPQ SC40 - DEC using standard

ifneq (,$(findstring alphaev6-dec-osf5.1,$(MACHINE)-$(VENDOR)-$(OS)))
  PPFC            := f90
  FC            := f90
  PFC           := f90
  FFLAGS1       := $(INCDIRS) -O5 -fast -pipeline -align dcommons -assume byterecl
  FFLAGS2       := $(FFLAGS1)
  FFLAGS3       := $(FFLAGS1)
  DA            := -DREAL8 -DCSCA
  DP            := -DREAL8 -DCMPI -DIBM -DCSCA
  DPRE          := -DREAL8 -DIBM
  IMODS         := -I
  CC            := cc
  CCBE          :=  $(CC)
  CFLAGS        := $(INCDIRS) -O2
  FLIBS          :=
  MSGLIBS       := -lmpi -lelan
  $(warning (INFO) Corresponding machine found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
endif

########################################################################
# Cray SV1 - Cray
# written by MEB 04/01/04 added by jgf45.06
ifneq (,$(findstring sv1-unicos,$(MACHINE)-$(OS)))
  PPFC           := f90
  FC            := f90
  PFC           := f90
  FFLAGS1	:=  $(INCDIRS) -dp -i32 -O2
  FFLAGS2	:=  $(INCDIRS) -dp -O2 -N 132
  FFLAGS3	:=  $(INCDIRS) -dp -O2 -N 132
  FIXED         :=  -f fixed
  FREE          :=  -f free
  DA  	        :=  -DREAL8 -DCRAY -DCSCA
  DP  	        :=  -DREAL8 -DCRAY -DCSCA -DCMPI
  DPRE	        :=  -DREAL8 -DCRAY
  IMODS		:=  -p
  CFLAGS	:=  $(INCDIRS) -I ../Lib -O2 -DCRAY
  FLIBS  	:=
  MSGLIBS	:=  -lmpi
  C_LDFLAGS     :=
  CCBE          :=  $(CC)
  $(warning (INFO) Corresponding machine found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
endif

########################################################################
# Cray-X1- Cray-X1
# written by MEB 04/01/04 added by jgf45.06
ifneq (,$(findstring x1-unicos,$(MACHINE)-$(OS)))
  PPFC          :=  ftn
  FC            :=  ftn
  PFC           :=  ftn
  FFLAGS1	:=  $(INCDIRS) -O2 -Oaggress
  FFLAGS2	:=  $(INCDIRS) -O2 -Oaggress
  FFLAGS3       :=  $(INCDIRS) -O2 -Oaggress
  FIXED         :=  -f fixed
  FREE          :=  -f free
  DA  	        :=  -DREAL8 -DCRAYX1 -DCVEC
  DP  	        :=  -DREAL8 -DCRAYX1 -DCVEC -DCMPI
  DPRE	        :=  -DREAL8 -DCRAYX1 -UCRAY
  IMODS		:=  -p
  CC            :=  cc
  CCBE          :=  $(CC)
  CFLAGS	:=  $(INCDIRS) -I ../Lib -O2 -DCRAYX1 -UCRAY
  FLIBS  	:=
  MSGLIBS	:=  -lmpi
  C_LDFLAGS     :=
  $(warning (INFO) Corresponding machine found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
endif

########################################################################
# SGI Origin
ifneq (,$(findstring mips-irix,$(MACHINE)-$(OS)))
  PPFC            := f90
  FC              := f90
  PFC             := f90
  FFLAGS1	  := $(INCDIRS) -O2 -OPT:Olimit=4257 -OPT:reorg_common=ON
  FFLAGS2	  := $(INCDIRS) -O2 -OPT:Olimit=4257 -OPT:reorg_common=ON
  FFLAGS3	  := $(INCDIRS) -O2 -OPT:Olimit=4257 -OPT:reorg_common=ON
  DA	          :=  -DREAL8 -DSGI -DCSCA
  DP	          :=  -DREAL8 -DSGI -DCSCA -DCMPI
  DPRE	          :=  -DREAL8 -DSGI
  IMODS   	  :=  -I
  CC              :=  cc
  CCBE            :=  $(CC)
  CFLAGS          :=  $(INCDIRS) -O2 -DSGI
  FLIBS		  :=
  MSGLIBS	  := -lmpi
  $(warning (INFO) Corresponding machine found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
endif
########################################################################


########################################################################
# powerpc-apple-darwin using absoft

ifneq (,$(findstring powerpc-darwin,$(MACHINE)-$(OS)))
  PPFC	        := f90
  FC	        := f90
  PFC	        := mpif77
  FFLAGS1	:=  $(INCDIRS) -w -O3 -m64 -cpu:g5 -f fixed -W132 -I . -DLINUX
  FFLAGS2	:=  $(INCDIRS) -w -O3 -m64 -cpu:g5 -N11 -f fixed -W132 -I .
  FFLAGS3	:=  $(INCDIRS) -w -O3 -m64 -cpu:g5 -N11 -f fixed -W132 -I .
  DA  	   	:=  -DREAL8 -DCSCA -DLINUX
  DP  	   	:=  -DREAL8 -DCSCA -DCMPI -DLINUX
  DPRE	   	:=  -DREAL8 -DLINUX
  IMODS  	:=  -p
  CC            :=  gcc
  CCBE          :=  $(CC)
  CFLAGS        :=  $(INCDIRS) -m64 -mpowerpc64 -O2 -DLINUX
  LDFLAGS	:=
  FLIBS	        :=  -lU77
  MSGLIBS	:=  -lm
  $(warning (INFO) Corresponding machine found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
endif
########################################################################


########################################################################
# i386-apple-darwin using intel

ifneq (,$(findstring i386-darwin,$(MACHINE)-$(OS)))
  PPFC	        := ifort
  FC	        := ifort
  PFC	        := mpif77
  FFLAGS1       :=  $(INCDIRS) -nowarn -O3    -fixed -132 -check all -traceback -DLINUX -DNETCDF_DEBUG -I .
# FFLAGS1	:=  $(INCDIRS) -nowarn -O3    -fixed -132 -DIBM -I .
  FFLAGS2	:=  $(INCDIRS) -nowarn -O3    -fixed -132 -I .
  FFLAGS3	:=  $(INCDIRS) -nowarn -O3    -fixed -132 -I .
  DA  	   	:=  -DREAL8 -DCSCA -DLINUX
  DP  	   	:=  -DREAL8 -DCSCA -DLINUX -DCMPI -DNETCDF_DEBUG
  DPRE	   	:=  -DREAL8 -DLINUX
  IMODS  	:=  -I
  CC            :=  gcc
  CCBE          :=  $(CC)
  CFLAGS        :=  $(INCDIRS) -O3 -DLINUX
  LDFLAGS	:=
  FLIBS	        :=
  MSGLIBS	:=
  $(warning (INFO) Corresponding machine found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
endif
########################################################################
ifneq ($(FOUND), TRUE)
     $(warning (WARNING) None of the platforms in cmplrflags.mk match your platform. As a result, the specific compilers and flags that are appropriate for you could not be specified. Please edit the cmplrflags.mk file to include your machine and operating system. Continuing with generic selections for compilers.)
  PPFC	        := f90
  FC	        := f90
  PFC	        := mpif90
  FFLAGS1	:=  $(INCDIRS)
  FFLAGS2	:=  $(FFLAGS1)
  FFLAGS3	:=  $(FFLAGS1)
  DA  	   	:=  -DREAL8 -DCSCA -DLINUX
  DP  	   	:=  -DREAL8 -DCSCA -DLINUX -DCMPI
  DPRE	   	:=  -DREAL8 -DLINUX
  IMODS  	:=  -I
  CC            :=  cc
  CCBE          :=  $(CC)
  CFLAGS        :=  $(INCDIRS) -DLINUX
  LDFLAGS	:=
  FLIBS	        :=
  MSGLIBS	:=
endif
ifeq ($(MULTIPLE),TRUE)
     $(warning (WARNING) More than one match in cmplrflags.mk. This may result in the wrong compilers being selected. Please check the cmplrflags.mk file to ensure that only one set of compiler flags is specified for your platform.)
endif
