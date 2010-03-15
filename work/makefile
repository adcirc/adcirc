SHELL:=/bin/sh
SRCDIR?=..
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

########################################################################
#  Get Canonical Machine NAME from config.guess
#
$(warning (INFO) Guessing the type of platform ADCIRC will run on...)

NAME     := $(shell ./config.guess)
LIST     := $(subst -, ,$(NAME))
MACHINE  := $(word 1, $(LIST))
VENDOR   := $(word 2, $(LIST))
OS       := $(subst  $(MACHINE)-$(VENDOR)-,,$(strip $(NAME)))

VTAG     := $(shell echo $(TAG) | sed -e 's/_/./g' -e 's/V//')

PARALLEL_ADCIRC   := padcirc
$(warning (INFO) Name is $(NAME), Machine is $(MACHINE), Vendor is $(VENDOR), and OS is $(OS).)
include cmplrflags.mk

ifneq ($(compiler),)
$(warning (INFO) The compiler variable in cmplrflags.mk is set to $(compiler).)
endif
$(warning (INFO) The following compilers have been selected...)
$(warning (INFO) The Fortran compiler for adcprep is set to $(PPFC).)
$(warning (INFO) The serial Fortran compiler is set to $(FC).)
$(warning (INFO) The parallel Fortran compiler is set to $(PFC).)
$(warning (INFO) The C compiler is set to $(CC).)

# Specify NETCDF=enable to gmake on the command line to activate netcdf i/o
ifeq ($(NETCDF),enable)
  PPFC    := $(PPFC) -DNETCDF -I$(NETCDFHOME)/include
  PFC     := $(PFC) -DNETCDF -I$(NETCDFHOME)/include
  FFLAGS1 := $(FFLAGS1) -DNETCDF -I$(NETCDFHOME)/include
  FFLAGS2 := $(FFLAGS2) -DNETCDF -I$(NETCDFHOME)/include
  FFLAGS3 := $(FFLAGS3) -DNETCDF -I$(NETCDFHOME)/include
  FFLAGS4 := $(FFLAGS4) -DNETCDF -I$(NETCDFHOME)/include
endif

####################### Target Specific Rules ###################################

#                                      adcprep   
ifeq ($(BUILDTYPE),adcprep)
  CF:= $(PPFC)
  O_DIR:=odir1/
  FFLAGS:= $(FFLAGS1) $(DPRE) $(IMODS) $(O_DIR) 
  VPATH :=  $(SRCDIR)/prep:$(SRCDIR):$(SRCDIR)/src:$(SRCDIR)/wind 
  ifeq ($(NETCDF),enable)
    LIBS  := -Lodir_metis -lmetis -L$(NETCDFHOME)/lib -lnetcdf $(FLIBS)
  else
    LIBS  := -Lodir_metis -lmetis $(FLIBS)
  endif
  PMSG_OBJ:=
endif
#                                      adcprep_be
ifeq ($(BUILDTYPE),adcprep_be)
  CF:= $(FC)
  O_DIR:=odir2/
  FFLAGS:= $(FFLAGS1) $(DPRE) $(IMODS) $(O_DIR) 
  VPATH :=  $(SRCDIR)/prep:$(SRCDIR)/src:$(SRCDIR):$(SRCDIR)/wind 
  ifeq ($(NETCDF),enable)
    LIBS  := -Lodir_metis -lmetis -L$(NETCDFHOME)/lib -lnetcdf $(FLIBS) 
  else
    LIBS  := -Lodir_metis_be -lmetis_be $(FLIBS)
  endif
  PMSG_OBJ:=
endif
#                                      adcirc   
ifeq ($(BUILDTYPE),adcirc)
  CF:= $(FC)
  O_DIR:=odir3/
  FFLAGS:= $(FFLAGS2) $(DA) $(IMODS) $(O_DIR)
  VPATH:=  $(SRCDIR)/src:$(SRCDIR)/wind:$(SRCDIR):$(SRCDIR)/prep
  ifeq ($(NETCDF),enable)
    LIBS  := -L$(NETCDFHOME)/lib -lnetcdf $(FLIBS)
  endif
  MSG_MOBJ:= $(O_DIR)globalio.o
endif
#                                      $(PARALLEL_ADCIRC)   
ifeq ($(BUILDTYPE),$(PARALLEL_ADCIRC))
  CF:= $(PFC)
  O_DIR:=odir4/
  FFLAGS:= $(FFLAGS3) $(DP) $(IMODS) $(O_DIR)
  VPATH :=  $(SRCDIR)/src:$(SRCDIR)/wind:$(SRCDIR):$(SRCDIR)/prep 
  ifeq ($(NETCDF),enable)
    LIBS  := -L$(NETCDFHOME)/lib -lnetcdf $(FLIBS)
  endif
  MSG_MOBJ:= $(O_DIR)messenger.o $(O_DIR)writer.o
endif
#                                      p15
ifeq ($(BUILDTYPE),p15)
  CF:= $(PPFC)
  O_DIR:=odir5/
  FFLAGS:= $(FFLAGS2) $(IMODS) $(O_DIR)
  VPATH :=  $(SRCDIR)/wind
  MSG_MOBJ:=
endif
#                                      owi22
ifeq ($(BUILDTYPE),owi22)
  CF:= $(PPFC)
  O_DIR:=odir6/
  FFLAGS:= $(FFLAGS2) $(IMODS) $(O_DIR)
  VPATH :=  $(SRCDIR)/wind
  MSG_MOBJ:=
endif
#                                      build13
ifeq ($(BUILDTYPE),build13)
  CF:= $(PPFC)
  O_DIR:=odir7/
   FFLAGS:= $(FFLAGS2) $(IMODS) $(O_DIR)
  VPATH :=  $(SRCDIR)/util
  MSG_MOBJ:=
endif
#                                      build12
ifeq ($(BUILDTYPE),build12)
  CF:= $(PPFC)
  O_DIR:=odir8/
  FFLAGS:= $(FFLAGS2) $(IMODS) $(O_DIR)
  VPATH :=  $(SRCDIR)/util
  MSG_MOBJ:=
endif
#
ifeq ($(BUILDTYPE),buildstwave23)
  CF:= $(PPFC)
  O_DIR:=odir9/
  FFLAGS:= $(FFLAGS2) $(IMODS) $(O_DIR)
  VPATH :=  $(SRCDIR)/util
  MSG_MOBJ:=
endif
#                                      hot2asc      
ifeq ($(BUILDTYPE),hot2asc)
  CF:= $(PPFC)
  O_DIR:=odir10/
  FFLAGS:= $(FFLAGS2) $(IMODS) $(O_DIR) $(LIBS)
  VPATH :=  $(SRCDIR)/util
  MSG_MOBJ:=
endif
#                                      inflate
ifeq ($(BUILDTYPE),inflate)
  CF:= $(PPFC)
  O_DIR:=odir11/
  FFLAGS:= $(FFLAGS2) $(IMODS) $(O_DIR) $(LIBS)
  VPATH :=  $(SRCDIR)/util
  MSG_MOBJ:=
endif
#                                      hstime
ifeq ($(BUILDTYPE),hstime)
  CF:= $(PPFC)
  O_DIR:=odir12/
  FFLAGS:= $(FFLAGS2) $(IMODS) $(O_DIR)
  VPATH :=  $(SRCDIR)/util
  MSG_MOBJ:=
endif
#                                      adcpost (for 3D and harmonic)
ifeq ($(BUILDTYPE),adcpost)
  CF:= $(PPFC)
  O_DIR:=odir13/
  FFLAGS:= $(FFLAGS1) $(DPRE) $(IMODS) $(O_DIR) 
  VPATH :=  $(SRCDIR)/prep 
endif
#                                      init
ifeq ($(BUILDTYPE),init)
  CF:= $(FC)
  O_DIR:=odir14/
  FFLAGS:= $(FFLAGS2) $(DPRE) $(IMODS) $(O_DIR)
  VPATH :=  $(SRCDIR)/util:$(SRCDIR)/src:$(SRCDIR)
  LIBS  := -L$(NETCDFHOME)/lib -lnetcdf $(FLIBS)
endif
#                                      aswip
ifeq ($(BUILDTYPE),aswip)
  CF:= $(PPFC)
  O_DIR:=odir3/
  FFLAGS:= $(FFLAGS2) $(IMODS) $(O_DIR)
  VPATH :=  $(SRCDIR)/wind
  MSG_MOBJ:=
endif
#
$(O_DIR):
	mkdir -p $@

################# Module Source, Object, & Mod Files ######################

PREP_MSRC  =  presizes.F 
PGLO_MSRC  =  pre_global.F
ADC_MSRC   =  version.F sizes.F global.F global_3dvs.F 
POST_MSRC  =  post_global.F
HARM_MSRC  =  harm.F  
PREC_MSRC  =  precision.F
CONS_MSRC  =  constants.F
FITP_MSRC  =  fitpack.F
VORT_MSRC  =  vortex.F
LSQ_MSRC   =  lsq.F
WIND_MSRC  =  wind.F owiwind.F rs2.F
SOLV_MSRC  =  itpackv.F
NA_MSRC    =  nodalattr.F
NC_MSRC    =  netcdfio.F
GIO_MSRC   =  globalio.F

PREP_MOBJ:= $(patsubst %.F, $(O_DIR)%.o, $(PREP_MSRC) )
POST_MOBJ:= $(patsubst %.F, $(O_DIR)%.o, $(POST_MSRC) )
PGLO_MOBJ:= $(patsubst %.F, $(O_DIR)%.o, $(PGLO_MSRC) )
SOLV_MOBJ:= $(patsubst %.F, $(O_DIR)%.o, $(SOLV_MSRC) )
ADC_MOBJ := $(patsubst %.F, $(O_DIR)%.o, $(ADC_MSRC)  )
HARM_MOBJ := $(patsubst %.F, $(O_DIR)%.o, $(HARM_MSRC)  )
PREC_MOBJ := $(patsubst %.F, $(O_DIR)%.o, $(PREC_MSRC)  )
CONS_MOBJ := $(patsubst %.F, $(O_DIR)%.o, $(CONS_MSRC)  )
FITP_MOBJ := $(patsubst %.F, $(O_DIR)%.o, $(FITP_MSRC)  )
VORT_MOBJ := $(patsubst %.F, $(O_DIR)%.o, $(VORT_MSRC)  )
LSQ_MOBJ := $(patsubst %.F, $(O_DIR)%.o, $(LSQ_MSRC)  )
WIND_MOBJ := $(patsubst %.F, $(O_DIR)%.o, $(WIND_MSRC)  )
NA_MOBJ  := $(patsubst %.F, $(O_DIR)%.o, $(NA_MSRC)  )
NC_MOBJ  := $(patsubst %.F, $(O_DIR)%.o, $(NC_MSRC)  )
GIO_MOBJ  := $(patsubst %.F, $(O_DIR)%.o, $(GIO_MSRC)  ) 

############################# Source & Object Files ######################
#meb 04/20/2004 - added machdep.F onto POST_SRC
#kmd 04/14/2009 - added init_cond.F into ADC_SRC

METIS_SRC  =  metis.F
PREP_SRC   =  adcprep.F decomp.F read_global.F prep.F interp.F machdep.F
ADC_SRC    =  adcirc.F read_input.F cstart.F hstart.F init_cond.F timestep.F vsmy.F transport.F write_output.F massbal.F
POST_SRC   =  adcpost.F post.F compare.F diffmerge.F machdep.F
P15_SRC    =  p15.F
OWI_SRC    =  owi22.F
B13_SRC    =  build13.F
B12_SRC    =  build12.F
B23_SRC    =  buildstwave23.F
H2A_SRC    =  hot2asc.F
INF_SRC    =  inflate.F
HST_SRC    =  hstime.F
INI_SRC    =  init.F
ASW_SRC    =  aswip_1.0.3.F

METIS_OBJ:= $(patsubst %.F, $(O_DIR)%.o, $(METIS_SRC) )
PREP_OBJ:= $(patsubst %.F, $(O_DIR)%.o, $(PREP_SRC) ) $(O_DIR)mkdir.o
POST_OBJ:= $(patsubst %.F, $(O_DIR)%.o, $(POST_SRC) ) 
ADC_OBJ:= $(patsubst %.F, $(O_DIR)%.o, $(ADC_SRC) ) $(O_DIR)mkdir.o
P15_OBJ:= $(patsubst %.F, $(O_DIR)%.o, $(P15_SRC) )
OWI_OBJ:= $(patsubst %.F, $(O_DIR)%.o, $(OWI_SRC) )
B13_OBJ:= $(patsubst %.F, $(O_DIR)%.o, $(B13_SRC) )
B12_OBJ:= $(patsubst %.F, $(O_DIR)%.o, $(B12_SRC) )
B23_OBJ:= $(patsubst %.F, $(O_DIR)%.o, $(B23_SRC) )
H2A_OBJ:= $(patsubst %.F, $(O_DIR)%.o, $(H2A_SRC) )
INF_OBJ:= $(patsubst %.F, $(O_DIR)%.o, $(INF_SRC) )
HST_OBJ:= $(patsubst %.F, $(O_DIR)%.o, $(HST_SRC) )
ASW_OBJ:= $(patsubst %.F, $(O_DIR)%.o, $(ASW_SRC) )

######################## compilation rules ###############################

$(O_DIR)%.o  : %.c
	$(COMPILE.c) $< -o $@
$(O_DIR)%.o  : %.F
	$(CF) -c $(FFLAGS) -o $@  $<
	if [ "`echo *.mod`" != '*.mod' ]; then mv *.mod $(O_DIR); fi

########################## Executable Targets ############################

.PHONY: all metis metis_be adcprep adcprep_be adcpost adcirc $(PARALLEL_ADCIRC)  p15 owi22 build13 build12 buildstwave23 hot2asc inflate hstime aswip

all :  metis adcprep $(BACKEND_EXEC) adcpost adcirc $(PARALLEL_ADCIRC)

ifeq ($(MAKELEVEL),0)
   odir_metis:
	mkdir $@
   odir_metis_be:
	mkdir $@
   metis: odir_metis
	$(MAKE) -C $(SRCDIR)/metis/Lib/ CC="$(CC)" CFLAGS="$(CFLAGS)" O_DIR="$(CURDIR)/odir_metis/" LIBMETIS="$(CURDIR)/odir_metis/libmetis.a" ARFLAGS="$(ARFLAGS)"
   metis_be: odir_metis_be
	$(MAKE) -C $(SRCDIR)/metis/Lib/ CC="$(CCBE)" CFLAGS="$(CFLAGS)" O_DIR="$(CURDIR)/odir_metis_be/" LIBMETIS="$(CURDIR)/odir_metis_be/libmetis_be.a" ARFLAGS="$(ARFLAGS)"
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
   init:
	$(MAKE) BUILDTYPE=init $@
   aswip:
	$(MAKE) BUILDTYPE=aswip $@
else
   adcprep::  $(O_DIR)
   adcprep_be::  $(O_DIR)
   adcirc::   $(O_DIR)
   adcpost::  $(O_DIR)
   $(PARALLEL_ADCIRC)::  $(O_DIR)
   p15::      $(O_DIR)
   owi22::    $(O_DIR)
   build13::  $(O_DIR)
   build12::  $(O_DIR)
   buildstwave23:: $(O_DIR)
   hot2asc:: $(O_DIR)
   inflate:: $(O_DIR)
   hstime::   $(O_DIR)
   init::   $(O_DIR)
   aswip::  $(O_DIR)

ifeq ($(NETCDF),enable)
   adcprep ::  $(NC_MOBJ) $(METIS_OBJ) $(PREP_OBJ) 
	$(CF) $(FFLAGS) -o $@ $(wildcard $(O_DIR)*.o) $(LIBS) $(LDFLAGS)
else
   adcprep ::  $(METIS_OBJ) $(PREP_OBJ) 
	$(CF) $(FFLAGS) -o $@ $(wildcard $(O_DIR)*.o) $(LIBS) $(LDFLAGS)
endif
   adcprep_be ::  $(METIS_OBJ) $(PREP_OBJ) 
	$(CF) $(FFLAGS) -o $@ $(wildcard $(O_DIR)*.o) $(LIBS) $(LDFLAGS)

   adcpost ::  $(POST_OBJ)
	$(CF) $(FFLAGS) -o $@ $(wildcard $(O_DIR)*.o) $(LIBS) $(LDFLAGS)

ifeq ($(NETCDF),enable)
   adcirc ::  $(NC_MOBJ) $(ADC_OBJ) $(O_DIR)driver.o
	$(CF) $(FFLAGS) -o $@ $(wildcard $(O_DIR)*.o) $(LIBS) $(LDFLAGS)
else
   adcirc ::  $(ADC_OBJ) $(O_DIR)driver.o
	$(CF) $(FFLAGS) -o $@ $(wildcard $(O_DIR)*.o) $(LIBS) $(LDFLAGS)
endif
   $(PARALLEL_ADCIRC) ::  $(NA_MOBJ) $(MSG_MOBJ) $(ADC_OBJ) $(O_DIR)driver.o
	$(CF) $(FFLAGS) -o $@ $(wildcard $(O_DIR)*.o) $(LIBS) $(MSGLIBS)

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
	$(CF) $(FFLAGS) -o $@ $(wildcard $(O_DIR)*.o)

   aswip ::  $(ASW_OBJ)
	$(CF) $(FFLAGS) -o $@ $(PREC_MOBJ) $(CONS_MOBJ) $(VORT_MOBJ) $(LSQ_MOBJ) $(FITP_MOBJ) $(ASW_OBJ)

endif

########################## Misc Commands ####################################

clean:
	rm -f odir*/*.o  odir*/*.mod sizes.o
clobber:
	rm -r -f odir* 
	rm -f adcprep adcprep_be adcirc $(PARALLEL_ADCIRC) p15 owi22 build13 \
           build12 buildstwave23 hot2asc inflate hstime sizes.o adcpost
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


echo:
	@echo VTAG: $(VTAG)

cvstag:
	ifneq ($(TAG),)
	  echo "      module version"                                        >  ../version.F; \
	  echo '      character(80), parameter :: ADC_VERSION = "'$(VTAG)'"' >> ../version.F; \
	  echo '      end module'                                            >> ../version.F; \
	  cd ..;                                                                              \
	  cvs ci -m'moving to TAG_VERSION $(VTAG)' version.F;                                 \
	  cvs tag $(TAG_FLAGS) $(TAG)
	else
	  @echo "To cvs tag do: make cvstag TAG=V?_?"
	endif



########################## Defining the DAG  #################################

#  adcprep modules

$(O_DIR)version.o     :  version.F
$(O_DIR)presize.o     :  presize.F
$(O_DIR)pre_global.o  :  pre_global.F  $(PREP_MOBJ) $(O_DIR)version.o
#  adcprep

$(O_DIR)adcprep.o     :  adcprep.F  $(PGLO_MOBJ) $(PMSG_OBJ)
$(O_DIR)decomp.o      :  decomp.F   $(PGLO_MOBJ)
$(O_DIR)read_global.o :  read_global.F  $(PGLO_MOBJ)
ifeq ($(NETCDF),enable)
$(O_DIR)prep.o        :  prep.F   $(PGLO_MOBJ) $(NC_MOBJ)
else
$(O_DIR)prep.o        :  prep.F   $(PGLO_MOBJ)
endif
$(O_DIR)interp.o      :  interp.F  
$(O_DIR)machdep.o     :  machdep.F
$(O_DIR)metis.o       :  metis.F $(PGLO_MOBJ)
$(O_DIR)mkdir.o       :  mkdir.c cfi.h

#  adcpost modules

$(O_DIR)post_global.o :  post_global.F $(POST_MOBJ)

#  adcpost

$(O_DIR)adcpost.o     :  adcpost.F $(POST_MOBJ)
$(O_DIR)post.o        :  post.F  $(POST_MOBJ)
$(O_DIR)compare.o     :  compare.F
$(O_DIR)diffmerge.o   :  diffmerge.F

#  asymmetric hurricane wind model

$(O_DIR)precision.o   :  precision.F
$(O_DIR)constants.o   :  constants.F $(PREC_MOBJ)
$(O_DIR)lsq.o         :  lsq.F $(PREC_MOBJ)
$(O_DIR)vortex.o      :  vortex.F $(PREC_MOBJ) $(CONS_MOBJ) $(LSQ_MOBJ) $(FITP_MOBJ)
$(O_DIR)aswip.o      :  aswip_1.0.3.F $(PREC_MOBJ) $(CONS_MOBJ) $(LSQ_MOBJ) $(FITP_MOBJ) $(VORT_MOBJ)

#  adcirc & $(PARALLEL_ADCIRC) modules

$(O_DIR)version.o     :  version.F
$(O_DIR)sizes.o       :  sizes.F
$(O_DIR)global.o      :  global.F  $(O_DIR)sizes.o
$(O_DIR)nodalattr.o   :  nodalattr.F  $(ADC_MOBJ)
$(O_DIR)messenger.o   :  messenger.F  $(ADC_MOBJ) 
$(O_DIR)writer.o      :  writer.F  $(O_DIR)messenger.o
$(O_DIR)globalio.o    :  globalio.F  $(ADC_MOBJ)
$(O_DIR)harm.o        :  harm.F
$(O_DIR)wind.o        :  wind.F  $(VORT_MOBJ) $(PREC_MOBJ) $(CONS_MOBJ) $(FITP_MOBJ)
$(O_DIR)owiwind.o     :  owiwind.F  $(VORT_MOBJ) $(PREC_MOBJ) $(CONS_MOBJ)
$(O_DIR)rs2.o         :  rs2.F  $(VORT_MOBJ) $(PREC_MOBJ) $(CONS_MOBJ)
$(O_DIR)itpackv.o     :  itpackv.F    $(ADC_MOBJ)
$(O_DIR)netcdfio.o    :  netcdfio.F  $(ADC_MOBJ) $(NA_MOBJ) $(HARM_MOBJ) $(WIND_MOBJ)

#  adcirc & $(PARALLEL_ADCIRC)

$(O_DIR)read_input.o  : read_input.F 
$(O_DIR)transport.o   : transport.F  
$(O_DIR)cstart.o      : cstart.F    
$(O_DIR)hstart.o      : hstart.F  
$(O_DIR)init_cond.o   : init_cond.F 
$(O_DIR)timestep.o    : timestep.F
$(O_DIR)write_output.o: write_output.F  
$(O_DIR)massbal.o     : massbal.F
ifeq ($(NETCDF),enable)
$(O_DIR)adcirc.o      : adcirc.F $(ADC_MOBJ) $(HARM_MOBJ) $(WIND_MOBJ) $(SOLV_MOBJ) $(NA_MOBJ) $(GIO_MOBJ) $(NC_MOBJ)
else
$(O_DIR)adcirc.o      : adcirc.F $(ADC_MOBJ) $(HARM_MOBJ) $(WIND_MOBJ) $(SOLV_MOBJ) $(NA_MOBJ) $(GIO_MOBJ)
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

# netcdf output file initializer

$(O_DIR)init.o       : init.F $(ADC_MOBJ) $(NC_MOBJ)


