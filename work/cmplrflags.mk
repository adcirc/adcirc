ifeq ($(compiler),intel-ND)
  PPFC          :=  ifort
  FC            :=  ifort
  PFC           :=  mpiifort
  INCDIRS       :=  $(INCDIRS)
#  INCDIRS       :=  $(INCDIRS) -I$(WGRIB2HOME) -I$(DATETIMEHOME)/include
  FFLAGS1       :=  $(INCDIRS) -w -O3 -assume byterecl -132 -assume buffered_io
  ifeq ($(DEBUG),full)
     FFLAGS1    :=  $(INCDIRS) -g -O0 -traceback -debug -check all -FI -assume byterecl -132 -DEBUG -DALL_TRACE -DFULL_STACK -DFLUSH_MESSAGES
  endif
  FFLAGS2       :=  $(FFLAGS1)
  FFLAGS3       :=  $(FFLAGS1)
  DA            :=  -DREAL8 -DLINUX -DCSCA
  DP            :=  -DREAL8 -DLINUX -DCSCA -DCMPI #-DPOWELL
  DPRE          :=  -DREAL8 -DLINUX -DADCSWAN
  ifeq ($(SWAN),enable)
     DPRE       := $(DPRE) -DADCSWAN
  endif
  IMODS         :=  -I
  CC            := icc
  CCBE          := $(CC)
  CFLAGS        := $(INCDIRS) -O3 -m64 -mcmodel=medium -DLINUX
  FLIBS         :=
  ifeq ($(GLOBAL),enable)
     WGRIB2HOME    := $(SRCDIR)/lib/grib2/lib/
     DATETIMEHOME  := $(SRCDIR)/lib/datetime-fortran/build/
     FLIBS         := -lwgrib2_api -lwgrib2 -ljasper -L$(WGRIB2HOME) -ldatetime -L$(DATETIMEHOME)lib/
  endif
  ifeq ($(DEBUG),full)
     CFLAGS     := $(INCDIRS) -g -O0 -m64 -march=k8 -mcmodel=medium -DLINUX
  endif
  ifeq ($(NETCDF),enable)
     HDF5HOME=/apps/hdf5/1.10.4/intel/16.1.150/lib
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
  NETCDFHOME=/apps/netcdf/4.6.1/intel/16.1.150
endif
