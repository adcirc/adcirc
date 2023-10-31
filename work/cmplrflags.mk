# SRCDIR is set in makefile or on the compile line
INCDIRS := -I. -I$(SRCDIR)/prep

########################################################################
# Compiler flags for Linux operating system on 64bit x86 CPU
#
ifeq ($(MACHINE)-$(OS),x86_64-linux-gnu)
#
# ***NOTE*** User must select between various Linux setups
#            by commenting/uncommenting the appropriate compiler
#
#compiler=gnu
ifeq ($(compiler),gnu)
  PPFC		:=  gfortran
  FC		:=  gfortran
  PFC		:=  mpif90
  FFLAGS1	:=  $(INCDIRS) -O2 -mcmodel=medium -ffixed-line-length-none -march=corei7-avx -m64 -static
  FFLAGS2	:=  $(FFLAGS1)
  FFLAGS3	:=  $(FFLAGS1)
  DA		:=  -DREAL8 -DLINUX -DCSCA
  DP		:=  -DREAL8 -DLINUX -DCSCA -DCMPI -DHAVE_MPI_MOD
  DPRE		:=  -DREAL8 -DLINUX
  IMODS 	:=  -I
  CC		:= gcc
  CCBE		:= $(CC)
  CFLAGS	:= $(INCDIRS) -O2 -mcmodel=medium -DLINUX -march=corei7-avx -m64 -static-libgcc
  CLIBS	:=
  LIBS		:=
  MSGLIBS	:=
  ifeq ($(NETCDF),enable)
        FLIBS          := $(FLIBS) -L$(HDF5HOME) -lhdf5 -lhdf5_fortran
  endif
  $(warning (INFO) Corresponding compilers and flags found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
endif
ifeq ($(compiler),ampi)
  PPFC		:=  ampif90 -pieglobals -DAMPI_COMP -module CommonLBs
  FC		:=  ampif90 -pieglobals -DAMPI_COMP -module CommonLBs
  PFC		:=  ampif90 -pieglobals -DAMPI_COMP -module CommonLBs
#  for easier debugging, without migration, in p=vp mode with charm
#  built non-smp
#  PPFC		:=  ampif90 -DAMPI_COMP
#  FC		:=  ampif90 -DAMPI_COMP
#  PFC		:=  ampif90 -DAMPI_COMP
  FFLAGS1	:=  $(INCDIRS) -O2 -mcmodel=medium -ffixed-line-length-none -march=corei7-avx -m64

  ifeq ($(DEBUG),full)
    FFLAGS1	:=  $(INCDIRS) -g -O0 -ffixed-line-length-none -fbacktrace -fbounds-check -ffpe-trap=zero,invalid,overflow,denormal -DALL_TRACE -DFLUSH_MESSAGES -DFULL_STACK -DDEBUG_HOLLAND -DDEBUG_WARN_ELEV
  endif	
  FFLAGS2	:=  $(FFLAGS1)
  FFLAGS3	:=  $(FFLAGS1)
  DA		:=  -DREAL8 -DLINUX -DCSCA -DSKIP_DRY
  DP		:=  -DREAL8 -DLINUX -DCSCA -DCMPI -DHAVE_MPI_MOD -DSKIP_DRY #-DREPART
  DPRE		:=  -DREAL8 -DLINUX
  IMODS 	:=  -I
  CC		:= gcc
  CCBE		:= $(CC)
  CFLAGS	:= $(INCDIRS) -O2 -mcmodel=medium -DLINUX -march=corei7-avx -m64
  CLIBS	:=
  LIBS		:= 
  MSGLIBS	:=
  ifeq ($(NETCDF),enable)
        FLIBS          := $(FLIBS) -L$(HDF5HOME) -lhdf5 -lhdf5_fortran
  endif
  $(warning (INFO) Corresponding compilers and flags found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
endif


endif
#$(MACHINE)



########################################################################
# Compiler flags for Linux operating system on 64bit x86 CPU
#
ifeq ($(MACHINE)-$(OS),x86_64-linux-gnu)
#
# ***NOTE*** User must select between various Linux setups
#            by commenting/uncommenting the appropriate compiler
#
compiler=gnu
#compiler=g95
#compiler=gfortran
#compiler=intel
#compiler=intel-ND
#compiler=intel-lonestar
#compiler=intel-sgi
#compiler=cray_xt3
#compiler=cray_xt4
#compiler=cray_xt5
#compiler=pgi
#compiler=pgi-ranger
#compiler=diamond
#compiler=kraken
#compiler=utils
#compiler=xtintel
#compiler=circleci
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
  DA		:=  -DREAL8 -DLINUX -DCSCA -DSKIP_DRY
  DP		:=  -DREAL8 -DLINUX -DCSCA -DCMPI -DHAVE_MPI_MOD -DSKIP_DRY
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
    FFLAGS1	:=  $(INCDIRS) -g -O0 -ffixed-line-length-none -fbacktrace -fbounds-check -ffpe-trap=zero,invalid,overflow,denormal -DALL_TRACE -DFLUSH_MESSAGES -DFULL_STACK -DDEBUG_HOLLAND -DDEBUG_WARN_ELEV
  endif
  ifeq ($(DEBUG),compiler-warnings)
    FFLAGS1	:=  $(INCDIRS) -g -O0 -Wall -Wextra -ffixed-line-length-none -DALL_TRACE -DFLUSH_MESSAGES -DFULL_STACK
  endif
  ifeq ($(DEBUG),full-not-warnelev)
    FFLAGS1	:=  $(INCDIRS) -g -O0 -ffixed-line-length-none -fbacktrace -fbounds-check -ffpe-trap=zero,invalid,overflow,denormal -DALL_TRACE -DFLUSH_MESSAGES -DFULL_STACK -DDEBUG_HOLLAND
  endif
  ifeq ($(DEBUG),full-not-fpe)
    FFLAGS1	:=  $(INCDIRS) -g -O0 -ffixed-line-length-none -fbacktrace -fbounds-check -DALL_TRACE -DFLUSH_MESSAGES -DFULL_STACK -DDEBUG_HOLLAND
  endif
  ifeq ($(DEBUG),trace)
    FFLAGS1	:=  $(INCDIRS) -g -O0 -ffixed-line-length-none -fbacktrace -DALL_TRACE -DFLUSH_MESSAGES -DFULL_STACK
  endif
#  ifneq ($(MACHINENAME),jason-desktop)
#     FFLAGS1 := $(FFLAGS1) -fno-underscoring
#  endif
  FFLAGS2	:=  $(FFLAGS1)
  FFLAGS3	:=  $(FFLAGS1)
  DA		:=  -DREAL8 -DLINUX -DCSCA -DSKIP_DRY
  DP		:=  -DREAL8 -DLINUX -DCSCA -DCMPI -DSKIP_DRY
  DPRE		:=  -DREAL8 -DLINUX
  ifeq ($(SWAN),enable)
     DPRE               :=  -DREAL8 -DLINUX -DADCSWAN
  endif
  FLIBS         :=
  ifeq ($(NETCDF),enable)
     ifeq ($(MACHINENAME),penguin)
        # module purge
        # module load gcc/6.2.0 openmpi/2.1.2/gcc.6.2.0
        # curl -O ftp://ftp.unidata.ucar.edu/pub/netcdf/old/netcdf-4.2.1.1.tar.gz
        # curl -O ftp://ftp.unidata.ucar.edu/pub/netcdf/old/netcdf-fortran-4.2.tar.gz
        # curl -O https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.12/src/hdf5-1.8.12.tar.gz
        # CPPFLAGS=-I/home/jgflemin/local/include LDFLAGS=-L/home/jgflemin/local/lib ./configure --prefix=/home/jgflemin/local
        #
        # export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HOME}/local/lib
        # export PATH=${PATH}:${HOME}/local/bin
        NETCDFHOME := ${HOME}/local
        FFLAGS1 := $(FFLAGS1) -I${HOME}/local/include
        FFLAGS2 := $(FFLAGS2) -I${HOME}/local/include
        FFLAGS3 := $(FFLAGS3) -I${HOME}/local/include
        FLIBS := $(FLIBS) -L${HOME}/local/lib -lnetcdff -lnetcdf
     endif
     ifeq ($(MACHINENAME),jason-desktop)
        NETCDFHOME := /usr
        FFLAGS1 := $(FFLAGS1) -L/usr/lib/x86_64-linux-gnu
        FFLAGS2 := $(FFLAGS2) -L/usr/lib/x86_64-linux-gnu
        FFLAGS3 := $(FFLAGS3) -L/usr/lib/x86_64-linux-gnu
     endif
     ifeq ($(MACHINENAME),rostam)
        NETCDFHOME := /usr
        FFLAGS1 := $(FFLAGS1) -I/usr/lib64/gfortran/modules
        FFLAGS2 := $(FFLAGS2) -I/usr/lib64/gfortran/modules
        FFLAGS3 := $(FFLAGS3) -I/usr/lib64/gfortran/modules
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
  PFC           ?=  mpif90
  FFLAGS1       := $(INCDIRS) -O2 -g -traceback -FI -assume byterecl -132 -xSSE4.2 -assume buffered_io
  CFLAGS        := $(INCDIRS) -O2 -xSSE4.2 -m64 -mcmodel=medium -DLINUX
  FLIBS         :=
  ifeq ($(DEBUG),full)
     CFLAGS        := $(INCDIRS) -g -O0 -march=k8 -m64 -mcmodel=medium -DLINUX
  endif
  ifeq ($(DEBUG),full)
     FFLAGS1       :=  $(INCDIRS) -g -O0 -traceback -debug all -check all -ftrapuv -fpe0 -FI -assume byterecl -132 -DALL_TRACE -DFULL_STACK -DFLUSH_MESSAGES
  endif
  ifeq ($(DEBUG),full-not-fpe)
     FFLAGS1       :=  $(INCDIRS) -g -O0 -traceback -debug all -check all -FI -assume byterecl -132 -DALL_TRACE -DFULL_STACK -DFLUSH_MESSAGES
  endif
  ifeq ($(DEBUG),trace)
     FFLAGS1       :=  $(INCDIRS) -g -O0 -traceback -FI -assume byterecl -132 -DALL_TRACE -DFULL_STACK -DFLUSH_MESSAGES
  endif
  ifeq ($(DEBUG),buserror)
     FFLAGS1       :=  $(INCDIRS) -g -O0 -traceback -FI -assume byterecl -132 -DALL_TRACE -DFULL_STACK -DFLUSH_MESSAGES -check bounds
  endif
  ifeq ($(DEBUG),netcdf_trace)
     FFLAGS1       :=  $(INCDIRS) -g -O0 -traceback -FI -assume byterecl -132 -DNETCDF_TRACE -DFULL_STACK -DFLUSH_MESSAGES
  endif
  #
  ifeq ($(MACHINENAME),stampede2)
     FFLAGS1 := $(INCDIRS) -O3 -FI -assume byterecl -132 -xCORE-AVX2 -axCORE-AVX512,MIC-AVX512 -assume buffered_io
     CFLAGS  := $(INCDIRS) -O3 -DLINUX -xCORE-AVX2 -axCORE-AVX512,MIC-AVX512
     FLIBS   := $(INCDIRS) -xCORE-AVX2 -axCORE-AVX512,MIC-AVX512
     ifeq ($(DEBUG),trace)
        FFLAGS1 := $(INCDIRS) -g -O0 -traceback -FI -assume byterecl -132 -xCORE-AVX2 -axCORE-AVX512,MIC-AVX512 -assume buffered_io
        CFLAGS  := $(INCDIRS) -g -O0 -traceback -DLINUX -xCORE-AVX2 -axCORE-AVX512,MIC-AVX512
        FLIBS   := $(INCDIRS) -xCORE-AVX2 -axCORE-AVX512,MIC-AVX512
     endif
  endif
  ifeq ($(MACHINENAME),frontera) 
     FFLAGS1 := $(INCDIRS) -O3 -FI -assume byterecl -132 -xCORE-AVX512 -assume buffered_io
     CFLAGS  := $(INCDIRS) -O3 -DLINUX -xCORE-AVX512 
     FLIBS   := $(INCDIRS) -xCORE-AVX512 
     ifeq ($(DEBUG),trace)
        FFLAGS1 := $(INCDIRS) -g -O0 -traceback -FI -assume byterecl -132 -xCORE-AVX512 -assume buffered_io
        CFLAGS  := $(INCDIRS) -g -O0 -traceback -DLINUX -xCORE-AVX512
        FLIBS   := $(INCDIRS) -xCORE-AVX512 
     endif
  endif
  ifeq ($(MACHINENAME),queenbee)
     FFLAGS1 := $(INCDIRS) -O3 -FI -assume byterecl -132 -xSSE4.2 -assume buffered_io
     CFLAGS  := $(INCDIRS) -O3 -DLINUX -xSSE4.2
     FLIBS   := $(INCDIRS) -xSSE4.2
     ifeq ($(DEBUG),trace)
        FFLAGS1 := $(INCDIRS) -g -O0 -traceback -FI -assume byterecl -132 -xSSE4.2 -assume buffered_io
        CFLAGS  := $(INCDIRS) -g -O0 -traceback -DLINUX -xSSE4.2
        FLIBS   := $(INCDIRS) -xSSE4.2
     endif
  endif
  ifeq ($(MACHINENAME),supermic) 
     FFLAGS1 := $(INCDIRS) -O3 -FI -assume byterecl -132 -xAVX -assume buffered_io
     CFLAGS  := $(INCDIRS) -O3 -DLINUX -xAVX
     FLIBS   := $(INCDIRS) -xAVX
     ifeq ($(DEBUG),trace)
        FFLAGS1 := $(INCDIRS) -g -O0 -traceback -FI -assume byterecl -132 -xAVX -assume buffered_io
        CFLAGS  := $(INCDIRS) -g -O0 -traceback -DLINUX  -xAVX 
        FLIBS   := $(INCDIRS) -xAVX
     endif
  endif
  #
  #@jasonfleming Added to fix bus error on hatteras@renci
  ifeq ($(HEAP_ARRAYS),fix)
     FFLAGS1 := $(FFLAGS1) -heap-arrays unlimited
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
  CLIBS         :=
  MSGLIBS       :=
  ifeq ($(NETCDF),enable)
     ifeq ($(MACHINENAME),hatteras)
        NETCDFHOME  :=$(shell nf-config --prefix)
        FLIBS       :=$(FLIBS) -L$(NETCDFHOME)/lib -lnetcdff -lnetcdf
        FFLAGS1     :=$(FFLAGS1) -I$(NETCDFHOME)/include
        FFLAGS2     :=$(FFLAGS1)
        FFLAGS3     :=$(FFLAGS1)
     endif
     # jgf20150417 queenbee requires that the analyst load the netcdf and
     # netcdf_fortran modules prior to compiling or executing ADCIRC
     ifeq ($(MACHINENAME),queenbee)
        FLIBS       := $(FLIBS) -L/usr/local/packages/netcdf/4.2.1.1/INTEL-140-MVAPICH2-2.0/lib -lnetcdff -lnetcdf
        NETCDFHOME    :=/usr/local/packages/netcdf/4.2.1.1/INTEL-140-MVAPICH2-2.0
     endif
     ifeq ($(MACHINENAME),supermic)
        FLIBS      := $(FLIBS) -L/usr/local/packages/netcdf_fortran/4.2/INTEL-140-MVAPICH2-2.0/lib -lnetcdff -L/usr/local/packages/netcdf/4.2.1.1/INTEL-140-MVAPICH2-2.0/lib -lnetcdf -lnetcdf -liomp5 -lpthread
        NETCDFHOME :=/usr/local/packages/netcdf/4.2.1.1/INTEL-140-MVAPICH2-2.0/include
        FFLAGS1    :=$(FFLAGS1) -I/usr/local/packages/hdf5/1.8.12/INTEL-140-MVAPICH2-2.0/include
     endif
     ifeq ($(MACHINENAME),stampede)
        NETCDFHOME :=/opt/apps/intel17/netcdf/4.3.3.1/x86_64
        FLIBS      := $(FLIBS) -L$(NETCDFHOME)/lib -lnetcdff -lnetcdf
     endif
     ifeq ($(MACHINENAME),stampede2)
        NETCDFHOME :=/opt/apps/intel17/netcdf/4.3.3.1/x86_64
        FLIBS      := $(FLIBS) -L$(NETCDFHOME)/lib -lnetcdff -lnetcdf
        ifeq ($(USER),jgflemin)
           NETCDFHOME :=/work/00976/jgflemin/stampede2/local
           FLIBS      := $(FLIBS) -L$(NETCDFHOME)/lib -lnetcdff -lnetcdf
        endif
     endif
     ifeq ($(MACHINENAME),frontera)
        # specify NETCDFHOME on the command line or as an environment var
        FLIBS      := $(FLIBS) -L$(NETCDFHOME)/lib -lnetcdff -lnetcdf
     endif
     # @jasonfleming: Added support for lonestar5 at tacc.utexas.edu;
     # load the following module: netcdf/4.3.3.1
     ifeq ($(MACHINENAME),lonestar5)
        #NETCDFHOME :=/opt/apps/intel18/netcdf/4.3.3.1/x86_64
        # @jasonfleming: Updated support for lonestar5
        NETCDFHOME :=/opt/apps/intel18/netcdf/4.6.2/x86_64
        FLIBS      := $(FLIBS) -L$(NETCDFHOME)/lib -lnetcdff -lnetcdf
     endif
     # jgf20150817: Adding support for spirit.afrl.hpc.mil;
     # load the following modules: netcdf-fortran/intel/4.4.2
     # and hdf5/intel/1.8.12 and hdf5-mpi/intel/sgimpt/1.8.12
     ifeq ($(MACHINENAME),spirit)
        NETCDFHOME :=/app/wpostool/COST/netcdf-fortran-4.4.2/intel
        FLIBS      := $(FLIBS) -L/app/wpostool/COST/netcdf-c-4.3.1.1/intel/lib -L$(NETCDFHOME)/lib -L/app/wpostool/COST/hdf5-mpi/1.8.12/intel/sgimpt/lib -lnetcdff -lnetcdf
     endif
     # jgf20150420 mike requires that the analyst add netcdf to the softenv
     # with the following on the command line
     # soft add +netcdf-4.1.3-Intel-13.0.0
     ifeq ($(MACHINENAME),mike)
        FLIBS       := $(FLIBS) -L/usr/local/packages/netcdf/4.1.3/Intel-13.0.0/lib -lnetcdff -lnetcdf
        NETCDFHOME    :=/usr/local/packages/netcdf/4.1.3/Intel-13.0.0
     endif
     ifeq ($(MACHINENAME),killdevil)
        HDF5HOME       :=/nas02/apps/hdf5-1.8.5/lib
        NETCDFHOME     :=/nas02/apps/netcdf-4.1.1
        FLIBS          := $(FLIBS) -L$(HDF5HOME) -lhdf5 -lhdf5_fortran
     endif
  endif

# ------------
  ifeq ($(NETCDF),enable)
    ifeq ($(WDALTVAL),enable)
       DP  := $(DP) -DWDVAL_NETCDF 
    endif
  endif
# -----------

  #jgf20110519: For netcdf on topsail at UNC, use
  #NETCDFHOME=/ifs1/apps/netcdf/
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
  PPFC          :=  ifort
  FC            :=  ifort
  PFC           ?=  mpif90
#  FFLAGS1       :=  $(INCDIRS) -w -O3 -assume byterecl -132 -assume buffered_io #-i-dynamic

  ifeq ($(AMD),yes)
     FFLAGS1     :=  $(INCDIRS) -g -traceback -O2 -assume byterecl -132 -mcmodel=medium -shared-intel -assume buffered_io
  else
     FFLAGS1     :=  $(INCDIRS) -O3 -g -traceback -xSSE4.2 -assume byterecl -132 -mcmodel=medium -shared-intel -assume buffered_io
#      FFLAGS1     :=  $(INCDIRS) -O2 -g -traceback -xSSE4.2 -assume byterecl -132 -mcmodel=medium -shared-intel -assume buffered_io
#      FFLAGS1     :=  $(INCDIRS) -O0 -g -traceback -check bounds -xSSE4.2 -assume byterecl -132 -mcmodel=medium -shared-intel -assume buffered_io
  endif
  ifeq ($(DEBUG),full)
     FFLAGS1    :=  $(INCDIRS) -g -O0 -traceback -debug -check all -FI -assume byterecl -132 -DEBUG -DALL_TRACE -DFULL_STACK -DFLUSH_MESSAGES
  endif
  FFLAGS2       :=  $(FFLAGS1)
  FFLAGS3       :=  $(FFLAGS1)
  DA            :=  -DREAL8 -DLINUX -DCSCA -DSKIP_DRY
  DP            :=  -DREAL8 -DLINUX -DCSCA -DCMPI -DSKIP_DRY #-DNOFSBPG #-DNOIVB -DPOWELL
  DPRE          :=  -DREAL8 -DLINUX -DADCSWAN
  ifeq ($(SWAN),enable)
     DPRE       := $(DPRE) -DADCSWAN
  endif
  IMODS         :=  -I
  CC            := icc
  CCBE          := $(CC)
  CFLAGS        := $(INCDIRS) -O3 -m64 -mcmodel=medium -DLINUX
  FLIBS         := 
  ifeq ($(DATETIME),enable)
     DATETIMEHOME  := /pontus/gling/ADCIRC_Codes/adcirc-cg-GLOBAL_WP_20200513/lib/datetime-fortran/
     FLIBS         := -ldatetime -L$(DATETIMEHOME)lib/
  endif
  ifeq ($(GRIB2),enable)
     WGRIB2HOME    := /pontus/gling/ADCIRC_Codes/adcirc-cg-GLOBAL_WP_20200513/lib/grib2/lib/
     FLIBS         := $(FLIBS) -lwgrib2_api -lwgrib2 -ljasper -L$(WGRIB2HOME)
  endif
  ifeq ($(DEBUG),full)
     CFLAGS     := $(INCDIRS) -g -O0 -m64 -march=k8 -mcmodel=medium -DLINUX
  endif
  ifeq ($(NETCDF),enable)
     #HDF5HOME=/afs/crc.nd.edu/x86_64_linux/hdf/hdf5-1.8.6-linux-x86_64-static/lib
     HDF5HOME=/opt/crc/n/netcdf/4.7.0/intel/18.0      
     FLIBS      := $(FLIBS) -lnetcdff -L$(HDF5HOME) 
  endif
# ------------
  ifeq ($(NETCDF),enable)
    ifeq ($(WDALTVAL),enable)
       DP  := $(DP) -DWDVAL_NETCDF 
    endif
  endif
# -----------

  CLIBS         :=
  MSGLIBS       :=
  $(warning (INFO) Corresponding machine found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
  NETCDFHOME=/afs/crc.nd.edu/x86_64_linux/n/netcdf/4.7.0/intel/18.0/
endif

# SGI ICE X (e.g. topaz@ERDC) using Intel compilers, added by TCM
# jgf: Added flags for Thunder@AFRL.
ifeq ($(compiler),intel-sgi)
  PPFC          :=  ifort
  FC            :=  ifort
  PFC           ?=  mpif90
  CC            :=  icc -O2 -no-ipo
  CCBE          :=  icc -O2 -no-ipo
  FFLAGS1       :=  $(INCDIRS) -fixed -extend-source 132 -O2 -finline-limit=1000 -real-size 64 -no-ipo -assume buffered_io
#  FFLAGS1      :=  $(INCDIRS) -Mextend -g -O0 -traceback
  FFLAGS2       :=  $(FFLAGS1)
  FFLAGS3       :=  $(FFLAGS1) -assume buffered_stdout
  DA            :=  -DREAL8 -DLINUX -DCSCA
  DP            :=  -DREAL8 -DLINUX -DCMPI -DHAVE_MPI_MOD -DCSCA
  DPRE          :=  -DREAL8 -DLINUX
  CFLAGS        :=  $(INCDIRS) -DLINUX
  IMODS         :=  -module
  FLIBS         :=
  ifeq ($(NETCDF),enable)
     ifeq ($(MACHINENAME),topaz)
        NETCDFHOME  :=/apps/unsupported/netcdf/4.3.3.1-intel-15.0.3
     endif
     ifeq ($(MACHINENAME),thunder)
        # add the following lines to ~/.personal.bashrc:
        # module load costinit
        # module load git
        # module load netcdf-fortran/intel/4.4.2
        NETCDFHOME  :=/app/COST/netcdf-fortran/4.4.2/intel
     endif
     # for platforms other than topaz, specify NETCDFHOME on the command line
     FLIBS       := $(FLIBS) -lnetcdff
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
# Cray-XC30/40 (e.g. armstrong@NAVO& lightning@afrl) using Intel compilers, added by TCM
ifeq ($(compiler),xtintel)
  PPFC          :=  ftn
  FC            :=  ftn
  PFC           :=  ftn
  CC            :=  cc -O2 -no-ipo
  CCBE          :=  cc -O2 -no-ipo
  FFLAGS1       :=  $(INCDIRS) -fixed -extend-source 132 -O2 -default64 -finline-limit=1000 -real-size 64 -no-ipo -assume buffered_io
#  FFLAGS1      :=  $(INCDIRS) -Mextend -g -O0 -traceback
  FFLAGS2       :=  $(FFLAGS1)
  FFLAGS3       :=  $(FFLAGS1) -assume buffered_stdout
  DA            :=  -DREAL8 -DLINUX -DCSCA
  DP            :=  -DREAL8 -DLINUX -DCMPI -DHAVE_MPI_MOD -DCSCA
  DPRE          :=  -DREAL8 -DLINUX
  CFLAGS        :=  $(INCDIRS) -DLINUX
  IMODS         :=  -module
  FLIBS         :=
# When compiling with netCDF support, the HDF5 libraries must also
# be linked in, so the user must specify HDF5HOME on the command line.
# jgf20090518: on Jade, NETCDFHOME=/usr/local/usp/PETtools/CE/pkgs/netcdf-4.0
# jgf20090518: on Jade, HDF5HOME=${PET_HOME}/pkgs/hdf5-1.8.2/lib
  ifeq ($(NETCDF),enable)
     FLIBS          := $(FLIBS) -L$(HDF5HOME) -lhdf5 -lhdf5_fortran
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
#
#
# Compiler Flags for CircleCI Build Server
ifeq ($(compiler),circleci)
  PPFC		:=  ifort
  FC		:=  ifort
  PFC		:=  mpif90
  FFLAGS1	:=  $(INCDIRS) -O0 -132
  FFLAGS2	:=  $(FFLAGS1)
  FFLAGS3	:=  $(FFLAGS1)
  DA		:=  -DREAL8 -DLINUX -DCSCA
  DP		:=  -DREAL8 -DLINUX -DCSCA -DCMPI -DHAVE_MPI_MOD
  DPRE		:=  -DREAL8 -DLINUX -DADCSWAN
  IMODS 	:=  -I
  CC		:= icx
  CCBE		:= $(CC)
  CFLAGS	:= $(INCDIRS) -O0 -mcmodel=medium -DLINUX -m64
  CLIBS	:=
  LIBS		:=
  MSGLIBS	:=
  ifeq ($(NETCDF),enable)
     FLIBS          := $(FLIBS) -L${NETCDF_C_HOME}/lib -lnetcdff
  endif
  $(warning (INFO) Corresponding compilers and flags found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
endif
#
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
  PFC	        ?=  mpif90
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
    FFLAGS1	:=  $(INCDIRS) -g -O0 -ffixed-line-length-none -fbacktrace -fbounds-check -ffpe-trap=zero,invalid,overflow,denormal -DALL_TRACE -DFLUSH_MESSAGES -DFULL_STACK -DDEBUG_HOLLAND -DDEBUG_WARN_ELEV
  endif
  ifeq ($(DEBUG),compiler-warnings)
    FFLAGS1	:=  $(INCDIRS) -g -O0 -Wall -Wextra -ffixed-line-length-none -DALL_TRACE -DFLUSH_MESSAGES -DFULL_STACK
  endif
  ifeq ($(DEBUG),full-not-warnelev)
    FFLAGS1	:=  $(INCDIRS) -g -O0 -ffixed-line-length-none -fbacktrace -fbounds-check -ffpe-trap=zero,invalid,overflow,denormal -DALL_TRACE -DFLUSH_MESSAGES -DFULL_STACK -DDEBUG_HOLLAND
  endif
  ifeq ($(DEBUG),full-not-fpe)
    FFLAGS1	:=  $(INCDIRS) -g -O0 -ffixed-line-length-none -fbacktrace -fbounds-check -DALL_TRACE -DFLUSH_MESSAGES -DFULL_STACK -DDEBUG_HOLLAND
  endif
  ifeq ($(DEBUG),trace)
    FFLAGS1	:=  $(INCDIRS) -g -O0 -ffixed-line-length-none -fbacktrace -DALL_TRACE -DFLUSH_MESSAGES -DFULL_STACK
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
  FFLAGS1	:=  $(INCDIRS) -w -O3 -m64 -cpu:g5 -f fixed -W132 -I. -DLINUX
  FFLAGS2	:=  $(INCDIRS) -w -O3 -m64 -cpu:g5 -N11 -f fixed -W132 -I.
  FFLAGS3	:=  $(INCDIRS) -w -O3 -m64 -cpu:g5 -N11 -f fixed -W132 -I.
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
  FFLAGS1       :=  $(INCDIRS) -nowarn -O3    -fixed -132 -check all -traceback -DLINUX -DNETCDF_DEBUG -I.
# FFLAGS1	:=  $(INCDIRS) -nowarn -O3    -fixed -132 -DIBM -I.
  FFLAGS2	:=  $(INCDIRS) -nowarn -O3    -fixed -132 -I.
  FFLAGS3	:=  $(INCDIRS) -nowarn -O3    -fixed -132 -I.
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
