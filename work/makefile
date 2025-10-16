SHELL:=/bin/bash
PWD=$(shell pwd)
SRCDIR?=$(shell dirname $(PWD))
#
#  Makefile to Build PADCIRC and its pre-processor
#  University of Texas's parallel version of the Hydrodynamics
#  Modeling program ADCIRC by J.J. Westerink and R.A. Luettich
#  vjp 3/29/2001
#  updated JJW 3/17/2003
#  updated for 3D JGF spring 2003
#  meb 4/20/2004
#  vjp 9/18/2006  made mods for Robert McLay's globalio module
#  vjp 9/19/2006  made mods for Shintaro's winds and waves
#  tcm 03/30/2011 added owi_ice module v49.64.01

########################################################################
#  Get Canonical Machine NAME from config.guess
#
$(warning (INFO) Guessing the type of platform ADCIRC will run on...)

NAME     := $(shell ./config.guess)
LIST     := $(subst -, ,$(NAME))
MACHINE  := $(word 1, $(LIST))
VENDOR   := $(word 2, $(LIST))
OS       := $(subst  $(MACHINE)-$(VENDOR)-,,$(strip $(NAME)))
PWD      := $(shell pwd)
VTAG     := $(shell echo $(TAG) | sed -e 's/_/./g' -e 's/V//')

PARALLEL_ADCIRC   := padcirc
LIBADC            := libadc.a
$(warning (INFO) Name is $(NAME), Machine is $(MACHINE), Vendor is $(VENDOR), and OS is $(OS).)
$(warning (INFO) The root directory for the build is $(SRCDIR))
$(warning (INFO) The ADCIRC version is $(shell ../scripts/adcirc_version.py --tag --create-version-file --directory $(SRCDIR)))
include cmplrflags.mk
#Casey 090302: Include variables/flags from SWAN.
include ../thirdparty/swan/macros.inc

ifneq ($(compiler),)
$(warning (INFO) The compiler variable in cmplrflags.mk is set to $(compiler).)
endif
$(warning (INFO) The following compilers have been selected...)
$(warning (INFO) The Fortran compiler for adcprep is set to $(PPFC).)
$(warning (INFO) The serial Fortran compiler is set to $(FC).)
$(warning (INFO) The parallel Fortran compiler is set to $(PFC).)
$(warning (INFO) The C compiler is set to $(CC).)
#ZC - Adding this as a warning to users to compile with the implicit array heap allocation options enabled
ifeq ($(FC),ifort)
ifneq ($(shell ulimit -s),unlimited)
    $(warning (WARNING) Intel compiler has been specified. The flag "--heap-arrays $(shell ulimit -s)" should be set)
    $(warning (WARNING) in both ADCIRC cmplrflags.mk FFLAGS1 and when compiling the netcdf-fortran library)
    $(warning (WARNING) using FCFLAGS to avoid potential issues with stack allocation of large implicit arrays.)
endif
endif

# Specify NETCDF=enable to gmake on the command line to activate netcdf i/o
ifeq ($(NETCDF),enable)
# jgf49.61: Don't add the -DADCNETCDF for the IBM compilers (see cmplrflags.mk)
ifeq ($(IBM),)
     PPFC    := $(PPFC) -DADCNETCDF -I$(NETCDFHOME)/include
     PFC     := $(PFC) -DADCNETCDF -I$(NETCDFHOME)/include
     FFLAGS1 := $(FFLAGS1) -DADCNETCDF -I$(NETCDFHOME)/include
     FFLAGS2 := $(FFLAGS2) -DADCNETCDF -I$(NETCDFHOME)/include
     FFLAGS3 := $(FFLAGS3) -DADCNETCDF -I$(NETCDFHOME)/include
     FFLAGS4 := $(FFLAGS4) -DADCNETCDF -I$(NETCDFHOME)/include
endif
endif
ifeq ($(NETCDF4),enable)
     PPFC    := $(PPFC) -DHAVE_NETCDF4
     PFC     := $(PFC) -DHAVE_NETCDF4
     FFLAGS1 := $(FFLAGS1) -DHAVE_NETCDF4
     FFLAGS2 := $(FFLAGS2) -DHAVE_NETCDF4
     FFLAGS3 := $(FFLAGS3) -DHAVE_NETCDF4
     FFLAGS4 := $(FFLAGS4) -DHAVE_NETCDF4
endif
ifeq ($(NETCDF4_COMPRESSION),enable)
     PPFC    := $(PPFC) -DNETCDF_CAN_DEFLATE
     PFC     := $(PFC) -DNETCDF_CAN_DEFLATE
     FFLAGS1 := $(FFLAGS1) -DNETCDF_CAN_DEFLATE
     FFLAGS2 := $(FFLAGS2) -DNETCDF_CAN_DEFLATE
     FFLAGS3 := $(FFLAGS3) -DNETCDF_CAN_DEFLATE
     FFLAGS4 := $(FFLAGS4) -DNETCDF_CAN_DEFLATE
endif
ifeq ($(XDMF),enable)
     $(warning (INFO) enabling xdmf lib ...)
     PPFC    := $(PPFC) -DADCXDMF -I$(XDMFHOME)/include -I$(SRCDIR)/src
     PFC     := $(PFC) -DADCXDMF -I$(XDMFHOME)/include -I$(SRCDIR)/src
     FFLAGS1 := $(FFLAGS1) -DADCXDMF -I$(XDMFHOME)/include -I$(SRCDIR)/src
     FFLAGS2 := $(FFLAGS2) -DADCXDMF -I$(XDMFHOME)/include -I$(SRCDIR)/src
     FFLAGS3 := $(FFLAGS3) -DADCXDMF -I$(XDMFHOME)/include -I$(SRCDIR)/src
     FFLAGS4 := $(FFLAGS4) -DADCXDMF -I$(XDMFHOME)/include -I$(SRCDIR)/src
     $(warning (INFO) XDMF library path is $(XDMFHOME).)
endif

FFLAGSGSW :=  
ifeq ($(GRIB2),enable)
     PPFC    := $(PPFC) -DGRIB2API
     PFC     := $(PFC) -DGRIB2API
     FFLAGS1 := $(FFLAGS1) -DGRIB2API
     FFLAGS2 := $(FFLAGS2) -DGRIB2API
     FFLAGS3 := $(FFLAGS3) -DGRIB2API
     FFLAGS4 := $(FFLAGS4) -DGRIB2API
     FFLAGSGSW := $(FFLAGSGSW) -I$(WGRIB2HOME)
endif
ifeq ($(DATETIME),enable) 
     PPFC    := $(PPFC) -DDATETIME
     PFC     := $(PFC) -DDATETIME
     FFLAGS1 := $(FFLAGS1) -DDATETIME
     FFLAGS2 := $(FFLAGS2) -DDATETIME
     FFLAGS3 := $(FFLAGS3) -DDATETIME
     FFLAGS4 := $(FFLAGS4) -DDATETIME
     FFLAGSGSW := $(FFLAGSGSW) -I$(DATETIMEHOME)include
endif
####################### Target Specific Rules ###################################

#                                      adcprep
ifeq ($(BUILDTYPE),adcprep)
  CF:= $(PPFC)
  O_DIR:=$(SRCDIR)/work/odir1/
  FFLAGS:= $(FFLAGS1) $(DPRE) $(IMODS) $(O_DIR) $(FFLAGSGSW)
  VPATH :=  $(SRCDIR)/prep:$(SRCDIR)/thirdparty/KDTREE2:$(SRCDIR):$(SRCDIR)/src:$(SRCDIR)/wind
  ifeq ($(XDMF),enable)
     FLIBS    :=  -L$(XDMFHOME)/lib -lXdmfCore -lXdmfUtils -lXdmf $(FLIBS)
  endif
  ifeq ($(NETCDF),enable)
    LIBS  := -Lodir_metis -lmetis -L$(NETCDFHOME)/lib -lnetcdf -lnetcdff $(FLIBS)
  else
    LIBS  := -Lodir_metis -lmetis $(FLIBS)
  endif
  PMSG_OBJ:=
endif
#                                      adcprep_be
ifeq ($(BUILDTYPE),adcprep_be)
  CF:= $(FC)
  O_DIR:=$(SRCDIR)/work/odir2/
  FFLAGS:= $(FFLAGS1) $(DPRE) $(IMODS) $(O_DIR)
  VPATH :=  $(SRCDIR)/prep:$(SRCDIR)/thirdparty/KDTREE2:$(SRCDIR)/src:$(SRCDIR):$(SRCDIR)/wind
  ifeq ($(XDMF),enable)
     FLIBS    :=  -L$(XDMFHOME)/lib -lXdmfCore -lXdmfUtils -lXdmf $(FLIBS)
  endif
  ifeq ($(NETCDF),enable)
    LIBS  := -Lodir_metis -lmetis -L$(NETCDFHOME)/lib -lnetcdf -lnetcdff $(FLIBS)
  else
    LIBS  := -Lodir_metis_be -lmetis_be $(FLIBS)
  endif
  PMSG_OBJ:=
endif
#                                      adcirc
ifeq ($(BUILDTYPE),adcirc)
  CF:= $(FC)
  O_DIR:=$(SRCDIR)/work/odir3/
  FFLAGS:= $(FFLAGS2) $(DA) $(IMODS) $(O_DIR) $(FFLAGSGSW)
  VPATH:=  $(SRCDIR)/src:$(SRCDIR)/thirdparty/KDTREE2:$(SRCDIR)/wind:$(SRCDIR):$(SRCDIR)/prep
  ifeq ($(XDMF),enable)
     FLIBS    :=  -L$(XDMFHOME)/lib -lXdmfCore -lXdmfUtils -lXdmf $(FLIBS)
  endif
  ifeq ($(NETCDF),enable)
    LIBS  := -L$(NETCDFHOME)/lib -lnetcdf -lnetcdff $(FLIBS)
  endif
  #MSG_MOBJ:= $(O_DIR)globalio.o
  MSG_MOBJ:=
endif
#Casey 090302: Add rules for the serial SWAN+ADCIRC.
#                                      adcswan
ifeq ($(BUILDTYPE),adcswan)
  CF:= $(FC)
  O_DIR:=$(SRCDIR)/work/odir33/
  FFLAGS:= $(FFLAGS2) -DCSWAN $(DA) $(IMODS) $(O_DIR) $(FFLAGSGSW)
  VPATH:=  $(SRCDIR)/src:$(SRCDIR)/thirdparty/KDTREE2:$(SRCDIR)/wind:$(SRCDIR):$(SRCDIR)/prep:$(SRCDIR)/thirdparty/swan
  ifeq ($(XDMF),enable)
     FLIBS    :=  -L$(XDMFHOME)/lib -lXdmfCore -lXdmfUtils -lXdmf $(FLIBS)
  endif
  ifeq ($(NETCDF),enable)
    LIBS  := -L$(NETCDFHOME)/lib -lnetcdf -lnetcdff $(FLIBS)
  endif
  MSG_MOBJ:= $(O_DIR)globalio.o
  ifneq ($(wildcard $(SRCDIR)/work/.swanlastcompile),)
    include $(SRCDIR)/work/.swanlastcompile
    ifneq ($(SWANLASTCOMPILE),serial)
        $(warning (WARNING) Removing prior parallel SWAN source files)
        $(shell rm -f ../thirdparty/swan/*.f90 ../thirdparty/swan/*.f ../thirdparty/swan/*.for)
    endif
  endif
  $(shell echo "SWANLASTCOMPILE := serial" > $(SRCDIR)/work/.swanlastcompile)
  $(shell perl ../thirdparty/swan/switch.pl $(swch) -adcirc ../thirdparty/swan/*.ftn ../thirdparty/swan/*.ftn90)
endif
#                                      $(PARALLEL_ADCIRC)
ifeq ($(BUILDTYPE),$(PARALLEL_ADCIRC))
  CF:= $(PFC)
  O_DIR:=$(SRCDIR)/work/odir4/
  FFLAGS:= $(FFLAGS3) $(DP) $(IMODS) $(O_DIR) $(FFLAGSGSW)
  VPATH :=  $(SRCDIR)/src:$(SRCDIR)/thirdparty/KDTREE2:$(SRCDIR)/wind:$(SRCDIR):$(SRCDIR)/prep
  ifeq ($(XDMF),enable)
     FLIBS    :=  -L$(XDMFHOME)/lib -lXdmfCore -lXdmfUtils -lXdmf $(FLIBS)
  endif
  ifeq ($(NETCDF),enable)
    LIBS  := -L$(NETCDFHOME)/lib -lnetcdf -lnetcdff $(FLIBS)
  endif
  #MSG_MOBJ:= $(O_DIR)messenger.o $(O_DIR)writer.o
  MSG_MOBJ:= $(O_DIR)messenger.o
endif
#Casey 090302: Add rules for the parallel SWAN+ADCIRC.
#                                      padcswan
ifeq ($(BUILDTYPE),padcswan)
  CF:= $(PFC)
  O_DIR:=$(SRCDIR)/work/odir44/
  FFLAGS:= $(FFLAGS3) -DCSWAN $(DP) $(IMODS) $(O_DIR) $(FFLAGSGSW)
  VPATH :=  $(SRCDIR)/src:$(SRCDIR)/thirdparty/KDTREE2:$(SRCDIR)/wind:$(SRCDIR):$(SRCDIR)/prep:$(SRCDIR)/thirdparty/swan
  ifeq ($(XDMF),enable)
     FLIBS    :=  -L$(XDMFHOME)/lib -lXdmfCore -lXdmfUtils -lXdmf $(FLIBS)
  endif
  ifeq ($(NETCDF),enable)
    LIBS  := -L$(NETCDFHOME)/lib -lnetcdf -lnetcdff $(FLIBS)
  endif
  #MSG_MOBJ:= $(O_DIR)messenger.o $(O_DIR)writer.o
  MSG_MOBJ:=  
  ifneq ($(wildcard $(SRCDIR)/work/.swanlastcompile),)
    include $(SRCDIR)/work/.swanlastcompile
    ifneq ($(SWANLASTCOMPILE),parallel)
        $(warning (WARNING) Removing prior serial SWAN source files)
        $(shell rm -f ../thirdparty/swan/*.f90 ../thirdparty/swan/*.f ../thirdparty/swan/*.for)
    endif
  endif
  $(shell echo "SWANLASTCOMPILE := parallel" > $(SRCDIR)/work/.swanlastcompile)
  $(shell perl ../thirdparty/swan/switch.pl $(swch) -pun -adcirc ../thirdparty/swan/*.ftn ../thirdparty/swan/*.ftn90)
endif
#                                    libadc.a
ifeq ($(BUILDTYPE),$(LIBADC))
  CF:= $(PFC)
  O_DIR:=$(SRCDIR)/work/odircp/
  FFLAGS:= $(FFLAGS3) $(DP) $(IMODS) $(O_DIR)
  VPATH := $(SRCDIR):$(SRCDIR)/src:$(SRCDIR)/thirdparty/KDTREE2:$(SRCDIR)/wind:$(SRCDIR):$(SRCDIR)/prep
  ifeq ($(XDMF),enable)
     FLIBS    :=  -L$(XDMFHOME)/lib -lXdmfCore -lXdmfUtils -lXdmf $(FLIBS)
  endif
  ifeq ($(NETCDF),enable)
    LIBS  := -L$(NETCDFHOME)/lib -lnetcdf -lnetcdff $(FLIBS)
  endif
  MSG_MOBJ:= $(O_DIR)messenger.o
endif
#                                      p15
ifeq ($(BUILDTYPE),p15)
  CF:= $(PPFC)
  O_DIR:=$(SRCDIR)/work/odir5/
  FFLAGS:= $(FFLAGS2) $(IMODS) $(O_DIR)
  VPATH :=  $(SRCDIR)/wind
  MSG_MOBJ:=
endif
#                                      owi22
ifeq ($(BUILDTYPE),owi22)
  CF:= $(PPFC)
  O_DIR:=$(SRCDIR)/work/odir6/
  FFLAGS:= $(FFLAGS2) $(IMODS) $(O_DIR)
  VPATH :=  $(SRCDIR)/wind
  MSG_MOBJ:=
endif
#                                      build13
ifeq ($(BUILDTYPE),build13)
  CF:= $(PPFC)
  O_DIR:=$(SRCDIR)/work/odir7/
   FFLAGS:= $(FFLAGS2) $(IMODS) $(O_DIR)
  VPATH :=  $(SRCDIR)/util
  MSG_MOBJ:=
endif
#                                      build12
ifeq ($(BUILDTYPE),build12)
  CF:= $(PPFC)
  O_DIR:=$(SRCDIR)/work/odir8/
  FFLAGS:= $(FFLAGS2) $(IMODS) $(O_DIR)
  VPATH :=  $(SRCDIR)/util
  MSG_MOBJ:=
endif
#
ifeq ($(BUILDTYPE),buildstwave23)
  CF:= $(PPFC)
  O_DIR:=$(SRCDIR)/work/odir9/
  FFLAGS:= $(FFLAGS2) $(IMODS) $(O_DIR)
  VPATH :=  $(SRCDIR)/util
  MSG_MOBJ:=
endif
#                                      hot2asc
ifeq ($(BUILDTYPE),hot2asc)
  CF:= $(PPFC)
  O_DIR:=$(SRCDIR)/work/odir10/
  FFLAGS:= $(FFLAGS2) $(IMODS) $(O_DIR) $(LIBS)
  VPATH :=  $(SRCDIR)/util
  MSG_MOBJ:=
endif
#                                      inflate
ifeq ($(BUILDTYPE),inflate)
  CF:= $(PPFC)
  O_DIR:=$(SRCDIR)/work/odir11/
  FFLAGS:= $(FFLAGS2) $(IMODS) $(O_DIR) $(LIBS)
  VPATH :=  $(SRCDIR)/util
  MSG_MOBJ:=
endif
#                                      hstime
ifeq ($(BUILDTYPE),hstime)
  CF:= $(PPFC)
  O_DIR:=$(SRCDIR)/work/odir12/
  FFLAGS:= $(FFLAGS2) $(IMODS) $(O_DIR) $(DA)
  ifeq ($(NETCDF),enable)
    LIBS  := -L$(NETCDFHOME)/lib -lnetcdf -lnetcdff $(FLIBS)
  endif
  VPATH :=  $(SRCDIR)/util
  MSG_MOBJ:=
endif
#                                      adcpost (for 3D and harmonic)
ifeq ($(BUILDTYPE),adcpost)
  CF:= $(PPFC)
  O_DIR:=$(SRCDIR)/work/odir13/
  FFLAGS:= $(FFLAGS1) $(DPRE) $(IMODS) $(O_DIR)
  VPATH :=  $(SRCDIR)/prep:$(SRCDIR)/src:$(SRCDIR)/thirdparty/KDTREE2:$(SRCDIR)
endif
#                                      aswip
ifeq ($(BUILDTYPE),aswip)
  CF:= $(PPFC)
  O_DIR:=$(SRCDIR)/work/odir14/
  FFLAGS:= $(FFLAGS2) $(IMODS) $(O_DIR) $(DA)
  VPATH := $(SRCDIR)/src:$(SRCDIR)/thirdparty/KDTREE2:$(SRCDIR)/wind:$(SRCDIR):$(SRCDIR)/prep
  ifeq ($(XDMF),enable)
     FLIBS    :=  -L$(XDMFHOME)/lib -lXdmfCore -lXdmfUtils -lXdmf $(FLIBS)
  endif
  ifeq ($(NETCDF),enable)
    LIBS  := -L$(NETCDFHOME)/lib -lnetcdf -lnetcdff $(FLIBS)
  endif
  MSG_MOBJ:= $(O_DIR)globalio.o
endif
#                                      adccmp
ifeq ($(BUILDTYPE),adccmp)
  CF:= $(PPFC)
  O_DIR:=$(SRCDIR)/work/odir17/
  FFLAGS:= $(FFLAGS2) $(IMODS) $(O_DIR) $(DA)
  VPATH :=  $(SRCDIR)/util
  MSG_MOBJ:=
endif
#
$(O_DIR):
	mkdir -p $@



################# Module Source, Object, & Mod Files ######################
HASH_MSRC  =  hashtable.F90
WD_MSRC    =  wetdry.F
MOM_MSRC   =  momentum.F
GWCE_MSRC  =  gwce.F
SUBD_MSRC  =  subdomain.F
SUBP_MSRC  =  subprep.F
CONT_MSRC  =  control.F
MESH_MSRC  =  mesh.F
BND_MSRC   =  boundaries.F
PREP_MSRC  =  presizes.F
PGLO_MSRC  =  pre_global.F
3D_MSRC    =  global_3dvs.F
ADC_MSRC   =  version.F sizes.F constants.F90 global.F
POST_MSRC  =  post_global.F
HARM_MSRC  =  harm.F
VORT_MSRC  =  vortex.F
VW1D_MSRC  =  vew1d.F90

ifeq ($(NETCDF),enable)
WIND_MSRC  =  wind.F owiwind.F rs2.F owi_ice.F owiwind_netcdf.F
else
WIND_MSRC  =  wind.F owiwind.F rs2.F owi_ice.F
endif

SOLV_MSRC  =  itpackv.F
NA_MSRC    =  nodalattr.F
SG_MSRC    =  subgridLookup.F
NC_MSRC    =  netcdfio.F90
NCE_MSRC   =  netcdf_error.F90
GIO_MSRC   =  globalio.F
WRIT_MSRC  =  write_output.F
WRITER_MSRC = writer.F 
NWS08_MSRC = nws08.F90
#Casey 090302: Add file for coupling to unstructured SWAN.
COUP_MSRC  =  couple2swan.F
#tcm v49.48.01 Adding File for kdtree2 fast search algorithm
KDTREE_MSRC = kdtree2.F
XDMF_MSRC = xdmfio.F
# DW sponge layer
SPG_MSRC = quadrature.F90 sponge_layer.F90
# WJP 3D baroclinic coupling
BC3D_MSRC = couple2baroclinic3D.F
# CPB baroclinic internal tide updates
GL2LOC_MSRC = gl2loc_mapping.F90
TIDE_MSRC = astronomic.F90 sun.F90 moon.F90 sun_moon_system.F90 ephemerides.F90 tidalpotential.F90
ITTIDE_MSRC = internaltide.F

HASH_MOBJ:= $(patsubst %.F90, $(O_DIR)%.o, $(HASH_MSRC) )
WD_MOBJ:= $(patsubst %.F, $(O_DIR)%.o, $(WD_MSRC) )
MOM_MOBJ:= $(patsubst %.F, $(O_DIR)%.o, $(MOM_MSRC) )
GWCE_MOBJ:= $(patsubst %.F, $(O_DIR)%.o, $(GWCE_MSRC) )
3D_MOBJ:= $(patsubst %.F, $(O_DIR)%.o, $(3D_MSRC) )
SUBD_MOBJ:= $(patsubst %.F, $(O_DIR)%.o, $(SUBD_MSRC) )
SUBP_MOBJ:= $(patsubst %.F, $(O_DIR)%.o, $(SUBP_MSRC) )
CONT_MOBJ:= $(patsubst %.F, $(O_DIR)%.o, $(CONT_MSRC) )
XDMF_MOBJ:= $(patsubst %.F, $(O_DIR)%.o, $(XDMF_MSRC) )
WRIT_MOBJ:= $(patsubst %.F, $(O_DIR)%.o, $(WRIT_MSRC) )
MESH_MOBJ:= $(patsubst %.F, $(O_DIR)%.o, $(MESH_MSRC) )
BND_MOBJ:= $(patsubst %.F, $(O_DIR)%.o, $(BND_MSRC) )
PREP_MOBJ:= $(patsubst %.F, $(O_DIR)%.o, $(PREP_MSRC) )
POST_MOBJ:= $(patsubst %.F, $(O_DIR)%.o, $(POST_MSRC) )
PGLO_MOBJ:= $(patsubst %.F, $(O_DIR)%.o, $(PGLO_MSRC) )
SOLV_MOBJ:= $(patsubst %.F, $(O_DIR)%.o, $(SOLV_MSRC) )
ADC_MOBJ := $(patsubst %.F, $(O_DIR)%.o, $(ADC_MSRC)  )
ADC_MOBJ := $(patsubst %.F90, $(O_DIR)%.o, $(ADC_MOBJ)  )
HARM_MOBJ := $(patsubst %.F, $(O_DIR)%.o, $(HARM_MSRC)  )
VORT_MOBJ := $(patsubst %.F, $(O_DIR)%.o, $(VORT_MSRC)  )
VW1D_MOBJ := $(patsubst %.F90, $(O_DIR)%.o, $(VW1D_MSRC)  )
WIND_MOBJ := $(patsubst %.F, $(O_DIR)%.o, $(WIND_MSRC)  )
NWS08_MOBJ := $(patsubst %.F90, $(O_DIR)%.o, $(NWS08_MSRC) )
NA_MOBJ  := $(patsubst %.F, $(O_DIR)%.o, $(NA_MSRC)  )
SG_MOBJ  := $(patsubst %.F, $(O_DIR)%.o, $(SG_MSRC)  )
NC_MOBJ  := $(patsubst %.F90, $(O_DIR)%.o, $(NC_MSRC)  )
GIO_MOBJ  := $(patsubst %.F, $(O_DIR)%.o, $(GIO_MSRC)  )
#Casey 090302: Add rules for coupling to unstructured SWAN.
COUP_MOBJ  := $(patsubst %.F, $(O_DIR)%.o, $(COUP_MSRC)  )
# tcm v49.48.01 Adding rule for kdtree2 fast search algorithm
KDTREE_MOBJ :=  $(patsubst %.F, $(O_DIR)%.o, $(KDTREE_MSRC) )
## DW
SPG_MOBJ :=  $(patsubst %.F90, $(O_DIR)%.o, $(SPG_MSRC) )
WRITER_MOBJ := $(patsubst %.F, $(O_DIR)%.o, $(WRITER_MSRC) )
# WJP
BC3D_MOBJ :=  $(patsubst %.F, $(O_DIR)%.o, $(BC3D_MSRC) )
# CPB
GL2LOC_MOBJ :=  $(patsubst %.F90, $(O_DIR)%.o, $(GL2LOC_MSRC) )
ITTIDE_MOBJ :=  $(patsubst %.F, $(O_DIR)%.o, $(ITTIDE_MSRC) )
TIDE_MOBJ := $(patsubst %.F90, $(O_DIR)%.o, $(TIDE_MSRC) )
NCE_MOBJ := $(patsubst %.F90, $(O_DIR)%.o, $(NCE_MSRC) )

############################# Source & Object Files ######################
#meb 04/20/2004 - added machdep.F onto POST_SRC

METIS_SRC  =  metis.F
PREP_SRC   =  adcprep.F decomp.F prep_weir.F read_global.F prep.F interp.F machdep.F
ADC_SRC    =  adcirc.F weir_boundary.F90 wetdry.F read_input.F normal_flow_boundary.F90 cstart.F hstart.F gwce.F momentum.F vew1d.F timestep.F vsmy.F transport.F
POST_SRC   =  adcpost.F post.F compare.F diffmerge.F machdep.F
P15_SRC    =  p15.F
OWI_SRC    =  owi22.F
B13_SRC    =  build13.F
B12_SRC    =  build12.F
B23_SRC    =  buildstwave23.F
H2A_SRC    =  hot2asc.F
INF_SRC    =  inflate.F
HST_SRC    =  hstime.F
ASW_SRC    =  aswip.F
CMP_SRC    =  adccmp.F

METIS_OBJ:= $(patsubst %.F, $(O_DIR)%.o, $(METIS_SRC) )
PREP_OBJ:= $(patsubst %.F, $(O_DIR)%.o, $(PREP_SRC) ) $(O_DIR)mkdir.o
POST_OBJ:= $(patsubst %.F, $(O_DIR)%.o, $(POST_SRC) )
ADC_OBJ:= $(patsubst %.F, $(O_DIR)%.o, $(ADC_SRC) ) $(O_DIR)mkdir.o
ADC_OBJ:= $(patsubst %.F90, $(O_DIR)%.o, $(ADC_OBJ) )
P15_OBJ:= $(patsubst %.F, $(O_DIR)%.o, $(P15_SRC) )
OWI_OBJ:= $(patsubst %.F, $(O_DIR)%.o, $(OWI_SRC) )
B13_OBJ:= $(patsubst %.F, $(O_DIR)%.o, $(B13_SRC) )
B12_OBJ:= $(patsubst %.F, $(O_DIR)%.o, $(B12_SRC) )
B23_OBJ:= $(patsubst %.F, $(O_DIR)%.o, $(B23_SRC) )
H2A_OBJ:= $(patsubst %.F, $(O_DIR)%.o, $(H2A_SRC) )
INF_OBJ:= $(patsubst %.F, $(O_DIR)%.o, $(INF_SRC) )
HST_OBJ:= $(patsubst %.F, $(O_DIR)%.o, $(HST_SRC) )
ASW_OBJ:= $(patsubst %.F, $(O_DIR)%.o, $(ASW_SRC) )
CMP_OBJ:= $(patsubst %.F, $(O_DIR)%.o, $(CMP_SRC) )

#################### SWAN Source & Object Files #########################
# Casey 090302: Added this section for coupling to unstructured SWAN.

SWAN_MSRC  = swmod1.f swmod2.f SwanSpectPart.f m_constants.f90 m_fileio.f90 \
serv_xnl4v5.f90 mod_xnl4v5.f90 \
SwanGriddata.f90 SwanGridobjects.f90 SwanCompdata.f90 SdsBabanin.f90 \
couple2adcirc.f90

SWAN_MOBJ := $(patsubst %.f,   $(O_DIR)%.o, $(SWAN_MSRC)  )
SWAN_MOBJ := $(patsubst %.f90, $(O_DIR)%.o, $(SWAN_MOBJ)  )

SWAN_SRC   = SwanQCM.f90 SwanIEM.f90 SwanBraggScat.f90 swanmain.f swanpre1.f fftpack51.f90 swanpre2.f \
swancom1.f swancom2.f swancom3.f swancom4.f \
swancom5.f swanout1.f swanout2.f swanser.f swanparll.f SwanReadGrid.f90 \
SwanReadADCGrid.f90 SwanReadTriangleGrid.f90 SwanReadEasymeshGrid.f90 \
SwanInitCompGrid.f90 SwanCheckGrid.f90 SwanCreateEdges.f90 SwanGridTopology.f90 SwanGridVert.f90 \
SwanGridCell.f90 SwanGridFace.f90 SwanPrintGridInfo.f90 SwanFindPoint.f90 \
SwanPointinMesh.f90 SwanBpntlist.f90 SwanPrepComp.f90 SwanVertlist.f90 SwanCompUnstruc.f90 \
SwanDispParm.f90 SwanPropvelX.f90 SwanSweepSel.f90 SwanPropvelS.f90 \
SwanTranspAc.f90 SwanTranspX.f90 SwanDiffPar.f90 SwanGSECorr.f90 SwanGradDepthorK.f90 \
SwanGradVel.f90 SwanInterpolatePoint.f90 \
SwanInterpolateAc.f90 SwanInterpolateOutput.f90 SwanConvAccur.f90 SwanConvStopc.f90 \
SwanFindObstacles.f90 SwanCrossObstacle.f90 \
SwanComputeForce.f90 SwanIntgratSpc.f90 SwanBndStruc.f90 SwanReadfort18.f90 \
SwanPunCollect.f90 \
SwanSumOverNodes.f90 SwanMinOverNodes.f90 SwanMaxOverNodes.f90 \
ocpids.f ocpcre.f ocpmix.f SwanVTKWriteHeader.f90 SwanVTKWriteData.f90 SwanVTKPDataSets.f90

SWAN_OBJ:= $(patsubst %.f,   $(O_DIR)%.o, $(SWAN_SRC) )
SWAN_OBJ:= $(patsubst %.f90, $(O_DIR)%.o, $(SWAN_OBJ) )

######################## compilation rules ###############################

$(O_DIR)%.o  : %.c
	$(COMPILE.c) $< -o $@
$(O_DIR)%.o  : %.F
	$(CF) -c $(FFLAGS) -o $@  $<
	if [ "`echo *.mod`" != '*.mod' ]; then mv *.mod $(O_DIR); fi
$(O_DIR)%.o  : %.F90
	$(CF) -c $(FFLAGS) -o $@  $<
	if [ "`echo *.mod`" != '*.mod' ]; then mv *.mod $(O_DIR); fi

# Casey 080106: Added the following rules for the unstructured SWAN.
$(O_DIR)%.o  : %.f
	$(CF) -c $(FLAGS_OPT) $(FLAGS_MSC) $(FLAGS_SER) $(IMODS) $(O_DIR) -o $@  $<
	if [ "`echo *.mod`" != '*.mod' ]; then mv *.mod $(O_DIR); fi
$(O_DIR)%.o  : %.f90
	$(CF) -c $(FLAGS_OPT) $(FLAGS90_MSC) $(FLAGS_SER) $(IMODS) $(O_DIR) -o $@  $<
	if [ "`echo *.mod`" != '*.mod' ]; then mv *.mod $(O_DIR); fi

########################## Executable Targets ############################

.PHONY: all libadc clean clobber metis metis_be adcprep adcprep_be adcpost adcirc $(PARALLEL_ADCIRC)  p15 owi22 build13 build12 buildstwave23 hot2asc inflate hstime aswip adccmp adcswan padcswan

all :  metis adcprep $(BACKEND_EXEC) adcpost adcirc $(PARALLEL_ADCIRC)

ifeq ($(MAKELEVEL),0)
   odir_metis:
	mkdir $@
   odir_metis_be:
	mkdir $@
   metis: odir_metis
	$(MAKE) -C $(SRCDIR)/thirdparty/metis/Lib/ CC="$(CC)" CFLAGS="$(CFLAGS)" O_DIR="$(CURDIR)/odir_metis/" LIBMETIS="$(CURDIR)/odir_metis/libmetis.a" ARFLAGS="$(ARFLAGS)"
   metis_be: odir_metis_be
	$(MAKE) -C $(SRCDIR)/thirdparty/metis/Lib/ CC="$(CCBE)" CFLAGS="$(CFLAGS)" O_DIR="$(CURDIR)/odir_metis_be/" LIBMETIS="$(CURDIR)/odir_metis_be/libmetis_be.a" ARFLAGS="$(ARFLAGS)"
   adcprep: metis
	$(MAKE) BUILDTYPE=adcprep  $@
   adcprep_be: metis_be
	$(MAKE) BUILDTYPE=adcprep_be CC="$(CCBE)" CFLAGS="$(CFLAGS)" $@
   adcpost:
	$(MAKE) BUILDTYPE=adcpost $@
   adcirc:
	$(MAKE) BUILDTYPE=adcirc CC="$(CCBE)" CFLAGS="$(CFLAGS)" $@
   $(PARALLEL_ADCIRC):
	$(MAKE) BUILDTYPE=$(PARALLEL_ADCIRC) CC="$(CCBE)" CFLAGS="$(CFLAGS)" $@
#Casey 090302: Added the following lines for coupling to unstructured SWAN.
   adcswan:
	$(MAKE) BUILDTYPE=adcswan CC="$(CCBE)" CFLAGS="$(CFLAGS)" $@
   padcswan:
	@perl ../thirdparty/swan/switch.pl $(swch) -pun -adcirc ../thirdparty/swan/*.ftn ../thirdparty/swan/*.ftn90
	$(MAKE) BUILDTYPE=padcswan CC="$(CCBE)" CFLAGS="$(CFLAGS)" $@
   libadc: $(LIBADC)
   $(LIBADC):
	@echo $(MAKE) LIBADC=$@ BUILDTYPE=$(LIBADC) CC="$(CCBE)" CFLAGS="$(CFLAGS)" $@
	$(MAKE) LIBADC=$@ BUILDTYPE=$(LIBADC) CC="$(CCBE)" CFLAGS="$(CFLAGS)" $@
   p15:
	$(MAKE) BUILDTYPE=p15 $@
   owi22:
	$(MAKE) BUILDTYPE=owi22 $@
   build13:
	$(MAKE) BUILDTYPE=build13 $@
   build12:
	$(MAKE) BUILDTYPE=build12 $@
   buildstwave23:
	$(MAKE) BUILDTYPE=buildstwave23 $@
   hot2asc:
	$(MAKE) BUILDTYPE=hot2asc $@
   inflate:
	$(MAKE) BUILDTYPE=inflate $@
   hstime:
	$(MAKE) BUILDTYPE=hstime $@
   aswip:
	$(MAKE) BUILDTYPE=aswip $@
   adccmp:
	$(MAKE) BUILDTYPE=adccmp $@
else
   adcprep::  $(O_DIR)
   adcprep_be::  $(O_DIR)
   adcirc::   $(O_DIR)
   adcpost::  $(O_DIR)
   $(PARALLEL_ADCIRC)::  $(O_DIR)
   $(LIBADC)::  $(O_DIR)
   p15::      $(O_DIR)
   owi22::    $(O_DIR)
   build13::  $(O_DIR)
   build12::  $(O_DIR)
   buildstwave23:: $(O_DIR)
   hot2asc:: $(O_DIR)
   inflate:: $(O_DIR)
   hstime::   $(O_DIR)
   aswip::  $(O_DIR)
   adccmp:  $(O_DIR)

ifeq ($(NETCDF),enable)
   adcprep ::  $(METIS_OBJ) $(NCE_MOBJ) $(KDTREE_MOBJ) $(SPG_MOBJ) $(BC3D_MOBJ) $(NC_MOBJ) $(WIND_MOBJ) $(NWS08_MOBJ) $(PREP_OBJ) $(MESH_MOBJ) $(NA_MOBJ) $(SUBP_MOBJ) $(HASH_MOBJ)
	$(CF) $(FFLAGS) -o $@ $(wildcard $(O_DIR)*.o) $(LIBS) $(LDFLAGS)
else
   adcprep ::  $(METIS_OBJ) $(NCE_MOBJ) $(KDTREE_MOBJ) $(SPG_MOBJ) $(BC3D_MOBJ) $(WIND_MOBJ) $(NWS08_MOBJ) $(PREP_OBJ) $(MESH_MOBJ) $(NA_MOBJ) $(SUBP_MOBJ) $(HASH_MOBJ)
	$(CF) $(FFLAGS) -o $@ $(wildcard $(O_DIR)*.o) $(LIBS) $(LDFLAGS)
endif
ifeq ($(NETCDF),enable)
   adcprep_be :: $(NC_MOBJ) $(METIS_OBJ) $(KDTREE_MOBJ) $(PREP_OBJ)
	$(CF) $(FFLAGS) -o $@ $(wildcard $(O_DIR)*.o) $(LIBS) $(LDFLAGS)
else
   adcprep_be ::  $(METIS_OBJ) $(KDTREE_MOBJ) $(PREP_OBJ)
	$(CF) $(FFLAGS) -o $@ $(wildcard $(O_DIR)*.o) $(LIBS) $(LDFLAGS)
endif
   adcpost ::  $(POST_OBJ) $(HARM_MOBJ)
	$(CF) $(FFLAGS) -o $@ $(wildcard $(O_DIR)*.o) $(LIBS) $(LDFLAGS)

ADCIRC_DEP := $(ADC_MOBJ) $(HAS_MOBJ) $(NCE_MOBJ) $(MESH_MOBJ) $(TIDE_MOBJ) $(ITTIDE_MOBJ) $(NA_MOBJ) $(KDTREE_MOBJ) $(SPG_MOBJ) $(BC3D_MOBJ) $(WIND_MOBJ) $(NWS08_MOBJ) $(ADC_OBJ) $(O_DIR)driver.o
ADCSWAN_DEP := $(ADC_MOBJ) $(HAS_MOBJ) $(NCE_MOBJ) $(MESH_MOBJ) $(TIDE_MOBJ) $(ITTIDE_MOBJ) $(NA_MOBJ) $(SWAN_MOBJ) $(KDTREE_MOBJ) $(SPG_MOBJ) $(BC3D_MOBJ) $(WIND_MOBJ) $(NWS08_MOBJ) $(ADC_OBJ) $(SWAN_OBJ) $(O_DIR)driver.o
ifeq ($(NETCDF),enable)
   ADCIRC_DEP := $(ADCIRC_DEP) $(NC_MOBJ)
   ADCSWAN_DEP := $(ADCSWAN_DEP) $(NC_MOBJ)
endif
ifeq ($(XDMF),enable)
   ADCIRC_DEP := $(ADCIRC_DEP) $(XDMF_MOBJ)
   ADCSWAN_DEP := $(ADCSWAN_DEP) $(XDMF_MOBJ)
endif

   adcirc ::  $(ADCIRC_DEP)
	$(CF) $(FFLAGS) -o $@ $(wildcard $(O_DIR)*.o) $(LIBS) $(LDFLAGS)

   $(PARALLEL_ADCIRC) :: $(ADC_MOBJ) $(HASH_MOBJ) $(MSG_MOBJ) $(NCE_MOBJ) $(MESH_MOBJ) $(TIDE_MOBJ) $(GL2LOC_MOBJ) $(ITTIDE_MOBJ) $(NA_MOBJ) $(KDTREE_MOBJ) $(SPG_MOBJ) $(BC3D_MOBJ) $(WIND_MOBJ) $(NWS08_MOBJ) $(ADC_OBJ) $(O_DIR)driver.o
	$(CF) $(FFLAGS) -o $@ $(wildcard $(O_DIR)*.o) $(LIBS) $(MSGLIBS)

ifeq ($(NETCDF),enable)
   $(LIBADC) :: $(O_DIR) $(ADC_MOBJ) $(HASH_MOBJ) $(MSG_MOBJ) $(NCE_MOBJ) $(MESH_MOBJ) $(TIDE_MOBJ) $(GL2LOC_MOBJ) $(ITTIDE_MOBJ) $(NA_MOBJ) \
	        $(KDTREE_MOBJ) $(SPG_MOBJ) $(BC3D_MOBJ) $(WIND_MOBJ) $(NWS08_MOBJ) $(GIO_MOBJ) $(AGIO_MOBJ) $(ADC_OBJ)
	ar $(ARFLAGS) $@  $(ADC_MOBJ) $(3D_MOBJ) $(MSG_MOBJ) $(HASH_MOBJ) $(TIDE_MOBJ) $(MESH_MOBJ) $(GL2LOC_MOBJ) $(ITTIDE_MOBJ) \
	         $(NA_MOBJ) $(GIO_MOBJ) $(AGIO_MOBJ) $(KDTREE_MOBJ) \
                 $(SPG_MOBJ) $(BC3D_MOBJ) $(WIND_MOBJ) $(NWS08_MOBJ) $(BND_MOBJ) $(HARM_MOBJ) $(WRIT_MOBJ) \
                 $(WRITER_MOBJ) $(SUBD_MOBJ) $(NC_MOBJ) $(VORT_MOBJ) $(SOLV_MOBJ) $(ADC_OBJ) 
	ar -ts $@
else
   $(LIBADC) :: $(O_DIR) $(ADC_MOBJ) $(HASH_MOBJ) $(MESH_MOBJ) $(TIDE_MOBJ) $(GL2LOC_MOBJ) $(ITTIDE_MOBJ) $(NA_MOBJ) $(MSG_MOBJ) \
	        $(KDTREE_MOBJ) $(SPG_MOBJ) $(BC3D_MOBJ) $(WIND_MOBJ) $(NWS08_MOBJ) $(GIO_MOBJ) $(AGIO_MOBJ) $(ADC_OBJ)
	ar $(ARFLAGS) $@  $(ADC_MOBJ) $(3D_MOBJ) $(MSG_MOBJ) $(HASH_MOBJ) $(MESH_MOBJ) $(TIDE_MOBJ) $(GL2LOC_MOBJ) $(ITTIDE_MOBJ) \
	         $(NA_MOBJ) $(GIO_MOBJ) $(AGIO_MOBJ) $(KDTREE_MOBJ) \
                 $(SPG_MOBJ) $(BC3D_MOBJ) $(WIND_MOBJ) $(NWS08_MOBJ) $(BND_MOBJ) $(HARM_MOBJ) $(WRIT_MOBJ) \
                 $(WRITER_MOBJ) $(SUBD_MOBJ) $(VORT_MOBJ) $(SOLV_MOBJ) $(ADC_OBJ) 
	ar -ts $@
endif

#Casey 080106: Added the following lines.
   adcswan :: $(O_DIR)
   adcswan :: $(ADCSWAN_DEP)
	$(CF) $(FFLAGS) -o $@ $(wildcard $(O_DIR)*.o) $(LIBS) $(LDFLAGS)
	if [ "`echo *.mod`" != '*.mod' ]; then mv *.mod $(O_DIR); fi

   padcswan :: $(O_DIR)
   padcswan :: $(NA_MOBJ) $(ADC_MOBJ) $(HASH_MOBJ) $(NCE_MOBJ) $(MESH_MOBJ) $(TIDE_MOBJ) $(GL2LOC_MOBJ) $(ITTIDE_MOBJ) $(MSG_MOBJ) $(SWAN_MOBJ) $(KDTREE_MOBJ) $(SPG_MOBJ) $(BC3D_MOBJ) $(WIND_MOBJ) $(NWS08_MOBJ) $(ADC_OBJ) $(SWAN_OBJ) $(O_DIR)driver.o
	$(CF) $(FFLAGS) -o $@ $(wildcard $(O_DIR)*.o) $(LIBS) $(MSGLIBS) $(LDFLAGS)
	if [ "`echo *.mod`" != '*.mod' ]; then mv *.mod $(O_DIR); fi

#ifeq ($(NETCDF),enable)
#   padcswan ::  $(NC_MOBJ) $(NA_MOBJ) $(MSG_MOBJ) $(SWAN_MOBJ) $(KDTREE_MOBJ) $(ADC_OBJ) $(SWAN_OBJ) $(O_DIR)driver.o
#	$(CF) $(FFLAGS) -o $@ $(wildcard $(O_DIR)*.o) $(LIBS) $(MSGLIBS) $(LDFLAGS)
#	if [ "`echo *.mod`" != '*.mod' ]; then mv *.mod $(O_DIR); fi
#else
#   padcswan ::  $(NA_MOBJ) $(MSG_MOBJ) $(SWAN_MOBJ) $(KDTREE_MOBJ) $(ADC_OBJ) $(SWAN_OBJ) $(O_DIR)driver.o
#	$(CF) $(FFLAGS) -o $@ $(wildcard $(O_DIR)*.o) $(LIBS) $(MSGLIBS) $(LDFLAGS)
#	if [ "`echo *.mod`" != '*.mod' ]; then mv *.mod $(O_DIR); fi
#endif

   p15     ::  $(P15_OBJ)
	$(CF) $(FFLAGS) -o $@ $(wildcard $(O_DIR)*.o)

   owi22   ::  $(OWI_OBJ)
	$(CF) $(FFLAGS) -o $@ $(wildcard $(O_DIR)*.o)

   build13 ::  $(B13_OBJ)
	$(CF) $(FFLAGS) -o $@ $(wildcard $(O_DIR)*.o)

   build12 ::  $(B12_OBJ)
	$(CF) $(FFLAGS) -o $@ $(wildcard $(O_DIR)*.o)

   buildstwave23 ::  $(B23_OBJ)
	$(CF) $(FFLAGS) -o $@ $(wildcard $(O_DIR)*.o)

   hot2asc ::  $(H2A_OBJ)
	$(CF) $(FFLAGS) -o $@ $(wildcard $(O_DIR)*.o)

   inflate ::  $(INF_OBJ)
	$(CF) $(FFLAGS) -o $@ $(wildcard $(O_DIR)*.o)

   hstime ::  $(HST_OBJ)
	$(CF) $(FFLAGS) -o $@ $(wildcard $(O_DIR)*.o) $(LIBS) $(LDFLAGS)

   aswip ::   $(ADCIRC_DEP) $(ASW_OBJ)
	$(CF) $(FFLAGS) -o $@ $(ADC_MOBJ) $(HAS_MOBJ) $(MESH_MOBJ) $(TIDE_MOBJ) $(ITTIDE_MOBJ) $(WIND_MOBJ) $(NWS08_MOBJ) $(NCE_MOBJ) $(NA_MOBJ) $(SG_MOBJ) $(PREC_MOBJ) $(CONS_MOBJ) $(VORT_MOBJ) $(LSQ_MOBJ) $(FITP_MOBJ) $(ASW_OBJ) $(KDTREE_MOBJ) $(HASH_MOBJ) $(BND_MOBJ) $(LIBS) $(LDFLAGS)

   adccmp :  $(O_DIR) $(CMP_OBJ)
	$(CF) $(FFLAGS) -o $@ $(wildcard $(O_DIR)*.o)


endif

########################## Misc Commands ####################################

clean:
	rm -f *.o odir*/*.o  *.mod odir*/*.mod actualflags.txt
	rm -f  ../thirdparty/swan/*.f ../thirdparty/swan/*.for ../thirdparty/swan/*.f90 .swanlastcompile
	rm -f ../thirdparty/nuopc/*.o ../thirdparty/nuopc/*.mod
	rm -f ../util/estofs_tide_fac/estofs_tide_fac

clobber: clean
	rm -r -f odir* util*
	rm -f adcprep adcprep_be adcirc $(PARALLEL_ADCIRC) adcswan padcswan p15 owi22 build13 aswip \
	      build12 buildstwave23 hot2asc inflate hstime constants.o sizes.o adcpost adccmp libadc.a

help:
	@echo "This makefile supports the following:"
	@echo "make all        - makes all seven targets"
	@echo "make adcprep    - makes the adcprep  front-end executable"
	@echo "make adcprep_be - makes the adcprep  back-end executable"
	@echo "make adcpost  - makes the adcpost  executable"
	@echo "make adcirc     - makes the serial adcirc executable"
	@echo "make $(PARALLEL_ADCIRC)  - makes the parallel adcirc executable"
	@echo "make p15        - makes the wind file generator for PBL"
	@echo "make owi22      - makes the wind file generator for OWI"
	@echo "make build13    - makes the fort.13 file generator"
	@echo "make build12    - makes the fort.12 file generator"
	@echo "make buildstwave23 - makes the fort.23 file generator"
	@echo "make hot2asc       - makes the hotstart debug program"
	@echo "make inflate       - makes the output file inflation program"
	@echo "make hstime   - makes the hot start time examiner"
	@echo "make adccmp     - makes the output file comparison program"
#Casey 091013: Add text for the SWAN+ADCIRC targets.
	@echo "make adcswan     - makes the serial SWAN+ADCIRC executable"
	@echo "make padcswan     - makes the parallel SWAN+ADCIRC executable"


########################## Defining the DAG  #################################

#  adcprep modules

$(O_DIR)version.o     :  version.F
$(O_DIR)presizes.o    :  presizes.F $(ADC_MOBJ) $(WIND_MOBJ) $(NWS08_MOBJ) $(O_DIR)global.o $(O_DIR)global_3dvs.o
$(O_DIR)kdtree2.o     :  kdtree2.F $(O_DIR)sizes.o
$(O_DIR)mesh.o        :  mesh.F $(O_DIR)version.o $(O_DIR)global.o $(O_DIR)boundaries.o $(O_DIR)hashtable.o
$(O_DIR)pre_global.o  :  pre_global.F  $(O_DIR)mesh.o $(KDTREE_MOBJ) $(PREP_MOBJ) $(O_DIR)global_3dvs.o
$(O_DIR)boundaries.o  :  boundaries.F $(O_DIR)sizes.o
$(O_DIR)subprep.o     :  subprep.F $(O_DIR)subprep.o $(O_DIR)global.o

#  adcprep

$(O_DIR)adcprep.o     :  adcprep.F  $(PGLO_MOBJ) $(PMSG_OBJ) $(ADC_MOBJ) $(SUBP_MOBJ)
$(O_DIR)decomp.o      :  decomp.F   $(PGLO_MOBJ)
$(O_DIR)read_global.o :  read_global.F  $(PGLO_MOBJ) $(KDTREE_MOBJ) $(ADC_MOBJ) $(NA_MOBJ) $(HARM_MOBJ) $(WIND_MOBJ) $(NWS08_MOBJ)
ifeq ($(NETCDF),enable)
$(O_DIR)prep.o        :  prep.F   $(PGLO_MOBJ) $(NC_MOBJ) $(HARM_MOBJ) $(NA_MOBJ)
else
$(O_DIR)prep.o        :  prep.F   $(PGLO_MOBJ) $(ADC_MOBJ) $(NA_MOBJ) $(HARM_MOBJ) $(WIND_MOBJ) $(NWS08_MOBJ) $(NA_MOBJ)
endif
$(O_DIR)interp.o      :  interp.F
$(O_DIR)machdep.o     :  machdep.F
$(O_DIR)metis.o       :  metis.F $(PGLO_MOBJ)
$(O_DIR)mkdir.o       :  mkdir.c cfi.h

#  adcpost modules

$(O_DIR)post_global.o :  post_global.F $(POST_MOBJ)

#  adcpost

$(O_DIR)adcpost.o     :  adcpost.F $(POST_MOBJ)
$(O_DIR)post.o        :  post.F  $(POST_MOBJ) $(HARM_MOBJ)
$(O_DIR)compare.o     :  compare.F
$(O_DIR)diffmerge.o   :  diffmerge.F

#  asymmetric hurricane wind model

$(O_DIR)vortex.o      :  vortex.F
$(O_DIR)aswip.o      :  aswip.F $(VORT_MOBJ) $(NA_MOBJ)

#  adcirc & $(PARALLEL_ADCIRC) modules

$(O_DIR)version.o     :  version.F
$(O_DIR)sizes.o       :  sizes.F
$(O_DIR)constants.o   :  constants.F90 $(O_DIR)sizes.o
$(O_DIR)global.o      :  global.F  $(O_DIR)sizes.o $(O_DIR)constants.o $(O_DIR)kdtree2.o
$(O_DIR)internaltide.o : internaltide.F $(O_DIR)sizes.o $(O_DIR)global.o
$(O_DIR)gl2loc_mapping.o : gl2loc_mapping.F90 $(O_DIR)sizes.o $(O_DIR)global.o

ifeq ($(BUILDTYPE),$(filter $(BUILDTYPE),$(PARALLEL_ADCIRC) padcswan))
$(O_DIR)nodalattr.o   :  nodalattr.F  $(O_DIR)sizes.o $(O_DIR)global.o $(O_DIR)messenger.o $(O_DIR)netcdf_error.o $(O_DIR)mesh.o $(O_DIR)hashtable.o $(O_DIR)messenger.o $(O_DIR)internaltide.o $(O_DIR)subgridLookup.o
else
$(O_DIR)nodalattr.o   :  nodalattr.F  $(O_DIR)sizes.o $(O_DIR)global.o $(O_DIR)netcdf_error.o $(O_DIR)mesh.o $(O_DIR)hashtable.o $(O_DIR)internaltide.o $(O_DIR)subgridLookup.o
endif

$(O_DIR)gwce.o        :  gwce.F $(O_DIR)itpackv.o $(O_DIR)nodalattr.o $(O_DIR)vew1d.o $(O_DIR)sponge_layer.o $(O_DIR)subdomain.o $(O_DIR)momentum.o
$(O_DIR)momentum.o    :  momentum.F $(O_DIR)vew1d.o $(O_DIR)subdomain.o
$(O_DIR)wetdry.o      :  wetdry.F $(O_DIR)subdomain.o
$(O_DIR)messenger.o   :  messenger.F90  $(O_DIR)sizes.o $(O_DIR)global.o $(O_DIR)global_3dvs.o
#!st3 100711: Added globalio.o to writer.F for HSWRITER MODULE
WRITER_DEP := writer.F  $(O_DIR)messenger.o $(O_DIR)globalio.o
ifeq ($(NETCDF),enable)
   WRITER_DEP := $(WRITER_DEP) $(O_DIR)netcdfio.o
endif
ifeq ($(XDMF),enable)
   WRITER_DEP := $(WRITER_DEP) $(O_DIR)xdmfio.o
endif
$(O_DIR)writer.o      : $(WRITER_DEP)
$(O_DIR)globalio.o    :  globalio.F  $(ADC_MOBJ) $(O_DIR)mesh.o $(O_DIR)global_3dvs.o
$(O_DIR)harm.o        :  harm.F $(O_DIR)mesh.o
$(O_DIR)wind.o        :  wind.F  $(VORT_MOBJ) $(MESH_MOBJ) $(HASH_MOBJ) $(NA_MOBJ) $(WIND_MOBJ) $(NWS08_MOBJ)
$(O_DIR)owiwind.o     :  owiwind.F  $(VORT_MOBJ)
$(O_DIR)owiwind_netcdf.o : owiwind_netcdf.F
$(O_DIR)owi_ice.o     :  owi_ice.F  $(VORT_MOBJ)
$(O_DIR)rs2.o         :  rs2.F  $(VORT_MOBJ)
ifeq ($(BUILDTYPE),$(filter $(BUILDTYPE),$(PARALLEL_ADCIRC) padcswan))
$(O_DIR)itpackv.o     :  itpackv.F $(O_DIR)sizes.o $(O_DIR)messenger.o
else
$(O_DIR)itpackv.o     :  itpackv.F $(O_DIR)sizes.o
endif
$(O_DIR)netcdfio.o    :  netcdfio.F90  $(ADC_MOBJ) $(NA_MOBJ) $(HARM_MOBJ) $(WIND_MOBJ) $(NWS08_MOBJ) $(GWCE_MOBJ) $(O_DIR)global_3dvs.o

#Casey 090302: Added the next line.
$(O_DIR)couple2swan.o :  couple2swan.F $(O_DIR)write_output.o
# jgf51.21.22: modularized mesh and boundaries
ifeq ($(BUILDTYPE),$(filter $(BUILDTYPE),$(LIBADC) $(PARALLEL_ADCIRC) padcswan))
$(O_DIR)mesh.o        :  mesh.F $(O_DIR)global.o $(O_DIR)boundaries.o $(O_DIR)messenger.o $(O_DIR)hashtable.o
else
$(O_DIR)mesh.o        :  mesh.F $(O_DIR)global.o $(O_DIR)boundaries.o $(O_DIR)hashtable.o
endif
$(O_DIR)boundaries.o  :  boundaries.F $(O_DIR)sizes.o
# jgf51.21.24: promoted write_output to a module

ifeq ($(BUILDTYPE),$(filter $(BUILDTYPE),$(LIBADC) $(PARALLEL_ADCIRC) padcswan))
WRITE_OUTPUT_DEP := write_output.F $(O_DIR)mesh.o $(O_DIR)global.o $(O_DIR)boundaries.o $(O_DIR)writer.o
else
WRITE_OUTPUT_DEP := write_output.F $(O_DIR)mesh.o $(O_DIR)global.o $(O_DIR)boundaries.o
endif
ifeq ($(NETCDF),enable)
   WRITE_OUTPUT_DEP := $(WRITE_OUTPUT_DEP) $(O_DIR)netcdfio.o
endif
ifeq ($(XDMF),enable)
   WRITE_OUTPUT_DEP := $(WRITE_OUTPUT_DEP) $(O_DIR)xdmfio.o
endif
$(O_DIR)write_output.o : $(WRITE_OUTPUT_DEP)
# jgf51.21.28: made xdmf support routines into a module
$(O_DIR)xdmfio.o      : xdmfio.F $(O_DIR)mesh.o $(O_DIR)global.o $(O_DIR)boundaries.o $(O_DIR)control.o
$(O_DIR)hashtable.o   : hashtable.F90
$(O_DIR)control.o     : control.F
$(O_DIR)subdomain.o   : subdomain.F $(O_DIR)mesh.o 
$(O_DIR)global_3dvs.o : global_3dvs.F $(O_DIR)global.o $(O_DIR)boundaries.o
# DW
$(O_DIR)sponge_layer.o : sponge_layer.F90 $(O_DIR)quadrature.o $(ADC_MOBJ)
$(O_DIR)quadrature.o : quadrature.F90
# WJP
$(O_DIR)couple2baroclinic3D.o : couple2baroclinic3D.F $(ADC_MOBJ) $(SPG_MOBJ)

#  adcirc & $(PARALLEL_ADCIRC)

$(O_DIR)read_input.o  : read_input.F $(SUBD_MOBJ)
$(O_DIR)transport.o   : transport.F
$(O_DIR)timestep.o    : timestep.F
ifeq ($(NETCDF),enable)
$(O_DIR)adcirc.o      : adcirc.F $(ADC_MOBJ) $(HARM_MOBJ) $(KDTREE_MOBJ) $(GWCE_MOBJ) $(MOM_MOBJ) $(WIND_MOBJ) $(NWS08_MOBJ) $(SOLV_MOBJ) $(NA_MOBJ) $(GIO_MOBJ) $(NC_MOBJ) $(COUP_MOBJ) $(WRIT_MOBJ)
else
$(O_DIR)adcirc.o      : adcirc.F $(ADC_MOBJ) $(HARM_MOBJ) $(KDTREE_MOBJ) $(GWCE_MOBJ) $(MOM_MOBJ) $(WIND_MOBJ) $(NWS08_MOBJ) $(SOLV_MOBJ) $(NA_MOBJ) $(GIO_MOBJ) $(COUP_MOBJ) $(WRIT_MOBJ)
endif
$(O_DIR)driver.o      : driver.F $(O_DIR)adcirc.o

# wind file generation

$(O_DIR)p15.o         : p15.F
$(O_DIR)owi22.o       : owi22.F

# fort.13 file generator

$(O_DIR)build13.o     : build13.F

# fort.12 file generator

$(O_DIR)build12.o     : build12.F

# fort.23 file generator

$(O_DIR)buildstwave23.o     : buildstwave23.F

# hotstart debug utility

$(O_DIR)hot2asc.o     : hot2asc.F

# converter from v46.32 output format to v47 format

$(O_DIR)inflate.o     : inflate.F

# hot start time examiner

$(O_DIR)hstime.o     : hstime.F

# output file comparison utility

$(O_DIR)adccmp.o      : adccmp.F

# Casey 090302: Added the following lines for SWAN.
$(O_DIR)swmod1.o                : swmod1.f
$(O_DIR)swmod2.o                : swmod2.f
$(O_DIR)SwanSpectPart.o         : SwanSpectPart.f
$(O_DIR)m_constants.o           : m_constants.f90
$(O_DIR)m_fileio.o              : m_fileio.f90
$(O_DIR)serv_xnl4v5.o           : serv_xnl4v5.f90
$(O_DIR)mod_xnl4v5.0            : mod_xnl4v5.f90
$(O_DIR)couple2adcirc.o         : couple2adcirc.f90
$(O_DIR)SwanGriddata.o          : SwanGriddata.f90
$(O_DIR)SwanGridobjects.o       : SwanGridobjects.f90
$(O_DIR)SwanCompdata.o          : SwanCompdata.f90
$(O_DIR)SdsBabanin.o            : SdsBabanin.f90
$(O_DIR)SwanIEM.o               : SwanIEM.f90
$(O_DIR)SwanQCM.o				: SwanQCM.f90
$(O_DIR)fftpack51.o             : fftpack51.f90
$(O_DIR)SwanBraggScat.o         : SwanBraggScat.f90
$(O_DIR)swanmain.o              : swanmain.f
$(O_DIR)swanpre1.o              : swanpre1.f
$(O_DIR)swanpre2.o              : swanpre2.f
$(O_DIR)swancom1.o              : swancom1.f
$(O_DIR)swancom2.o              : swancom2.f
$(O_DIR)swancom3.o              : swancom3.f
$(O_DIR)swancom4.o              : swancom4.f
$(O_DIR)swancom5.o              : swancom5.f
$(O_DIR)swanout1.o              : swanout1.f
$(O_DIR)swanout2.o              : swanout2.f
$(O_DIR)swanser.o               : swanser.f
$(O_DIR)swanparll.o             : swanparll.f
$(O_DIR)SwanReadGrid.o          : SwanReadGrid.f90
$(O_DIR)SwanReadADCGrid.o       : SwanReadADCGrid.f90
$(O_DIR)SwanReadTriangleGrid.o  : SwanReadTriangleGrid.f90
$(O_DIR)SwanReadEasymeshGrid.o  : SwanReadEasymeshGrid.f90
$(O_DIR)SwanInitCompGrid.o      : SwanInitCompGrid.f90
$(O_DIR)SwanCheckGrid.o         : SwanCheckGrid.f90
$(O_DIR)SwanCreateEdges.o       : SwanCreateEdges.f90
$(O_DIR)SwanGridTopology.o      : SwanGridTopology.f90
$(O_DIR)SwanGridVert.o          : SwanGridVert.f90
$(O_DIR)SwanGridCell.o          : SwanGridCell.f90
$(O_DIR)SwanGridFace.o          : SwanGridFace.f90
$(O_DIR)SwanPrintGridInfo.o     : SwanPrintGridInfo.f90
$(O_DIR)SwanFindPoint.o         : SwanFindPoint.f90
$(O_DIR)SwanPointinMesh.o       : SwanPointinMesh.f90
$(O_DIR)SwanBpntlist.o          : SwanBpntlist.f90
$(O_DIR)SwanPrepComp.o          : SwanPrepComp.f90
$(O_DIR)SwanVertlist.o          : SwanVertlist.f90
$(O_DIR)SwanCompUnstruc.o       : SwanCompUnstruc.f90
$(O_DIR)SwanDispParm.o          : SwanDispParm.f90
$(O_DIR)SwanPropvelX.o          : SwanPropvelX.f90
$(O_DIR)SwanSweepSel.o          : SwanSweepSel.f90
$(O_DIR)SwanPropvelS.o          : SwanPropvelS.f90
$(O_DIR)SwanTranspAc.o          : SwanTranspAc.f90
$(O_DIR)SwanTranspX.o           : SwanTranspX.f90
$(O_DIR)SwanDiffPar.o           : SwanDiffPar.f90
$(O_DIR)SwanGSECorr.o           : SwanGSECorr.f90
$(O_DIR)SwanGradDepthorK.o      : SwanGradDepthorK.f90
$(O_DIR)SwanGradVel.o           : SwanGradVel.f90
$(O_DIR)SwanInterpolatePoint.o  : SwanInterpolatePoint.f90
$(O_DIR)SwanInterpolateAc.o     : SwanInterpolateAc.f90
$(O_DIR)SwanInterpolateOutput.o : SwanInterpolateOutput.f90
$(O_DIR)SwanConvAccur.o         : SwanConvAccur.f90
$(O_DIR)SwanConvStopc.o         : SwanConvStopc.f90
$(O_DIR)SwanThreadBounds.o      : SwanThreadBounds.f90
$(O_DIR)SwanFindObstacles.o     : SwanFindObstacles.f90
$(O_DIR)SwanCrossObstacle.o     : SwanCrossObstacle.f90
$(O_DIR)SwanComputeForce.o      : SwanComputeForce.f90
$(O_DIR)SwanIntgratSpc.o        : SwanIntgratSpc.f90
$(O_DIR)SwanBndStruc.o          : SwanBndStruc.f90
$(O_DIR)SwanReadfort18.o        : SwanReadfort18.f90
$(O_DIR)SwanPunCollect.o        : SwanPunCollect.f90
$(O_DIR)SwanSumOverNodes.o      : SwanSumOverNodes.f90
$(O_DIR)SwanMinOverNodes.o      : SwanMinOverNodes.f90
$(O_DIR)SwanMaxOverNodes.o      : SwanMaxOverNodes.f90
$(O_DIR)ocpids.o                : ocpids.f
$(O_DIR)ocpcre.o                : ocpcre.f
$(O_DIR)ocpmix.o                : ocpmix.f

# jgf51.35: Write a file containing the actual compiler flags that were
# used so we can go back and verify later.
ifneq ($(BUILDTYPE),)
   FLAGSFILE := actualflags.txt
   $(shell echo ---------------------------------- >> $(FLAGSFILE))
   $(shell date >> $(FLAGSFILE))
   $(shell echo 'building $(BUILDTYPE)' >> $(FLAGSFILE))
   $(shell echo 'PPFC := $(PPFC)' >> $(FLAGSFILE))
   $(shell echo 'FC := $(FC)' >> $(FLAGSFILE))
   $(shell echo 'PFC := $(PFC)' >> $(FLAGSFILE))
   $(shell echo 'FFLAGS1 := $(FFLAGS1)' >> $(FLAGSFILE))
   $(shell echo 'FFLAGS2 := $(FFLAGS2)' >> $(FLAGSFILE))
   $(shell echo 'FFLAGS3 := $(FFLAGS3)' >> $(FLAGSFILE))
   $(shell echo 'DA := $(DA)' >> $(FLAGSFILE))
   $(shell echo 'DP := $(DP)' >> $(FLAGSFILE))
   $(shell echo 'DPRE := $(DPRE)' >> $(FLAGSFILE))
   $(shell echo 'IMODS := $(IMODS)' >> $(FLAGSFILE))
   $(shell echo 'CC := $(CC)' >> $(FLAGSFILE))
   $(shell echo 'CCBE := $(CCBE)' >> $(FLAGSFILE))
   $(shell echo 'CFLAGS := $(CFLAGS)' >> $(FLAGSFILE))
   $(shell echo 'CLIBS := $(CLIBS)' >> $(FLAGSFILE))
   $(shell echo 'LIBS := $(LIBS)' >> $(FLAGSFILE))
   $(shell echo 'MSGLIBS := $(MSGLIBS)' >> $(FLAGSFILE))
   # the following may have been supplied on the make command line
   $(shell echo 'compiler := $(compiler)' >> $(FLAGSFILE))
   $(shell echo 'MACHINENAME := $(MACHINENAME)' >> $(FLAGSFILE))
   $(shell echo 'NETCDF := $(NETCDF)' >> $(FLAGSFILE))
   $(shell echo 'NETCDFHOME := $(NETCDF)' >> $(FLAGSFILE))
   $(shell echo 'NETCDF4 := $(NETCDF4)' >> $(FLAGSFILE))
   $(shell echo 'NETCDF4_COMPRESSION := $(NETCDF4_COMPRESSION)' >> $(FLAGSFILE))
   $(shell echo 'XDMF := $(XDMFHOME)' >> $(FLAGSFILE))
   # get the compiler version while we're at it
   $(shell $(FC) --version | head -n 1 >> $(FLAGSFILE))
   $(shell echo ---------------------------------- >> $(FLAGSFILE))
endif
