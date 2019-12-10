################################################################################
##!! @author Saeed Moghimi (moghimis@gmail.com)
##!! @date 15/1/17 Original documentation

##! LOG !
##! 
################################################################################

ifneq ($(origin ESMFMKFILE), environment)
$(error Environment variable ESMFMKFILE was not set.)
endif

include $(ESMFMKFILE)




#ifneq ($(origin ADCDIR), environment)
#$(error Environment variable ADCDIR was not set.)
#endif

current_dir := $(shell pwd)
ADCDIR  := $(current_dir)/../..

INCDIRS := -I . -I$(ADCDIR)/work/odircp -I $(ADCDIR)/prep
LIBDIRS := -L . -L$(ADCDIR)/work 
LIBSS   := -ladc  -lnetcdf  -lnetcdff 

INCDIRS += -I${NETCDFHOME}/include/   
LIBDIRS += -L${NETCDFHOME}/lib/

INCDIRS += -I${HDF5HOME}/include/   
LIBDIRS += -L${HDF5HOME}/lib/

# unremark this line to define petlists
LIBDIRS += -DWITHPETLISTS_on 

# debug for gfortran
# LIBDIRS +=  -g -O0 -fbacktrace  -lhdf5 -lhdf5_fortran
################################################################################ $(INCDIRS) $(LIBDIRS)
################################################################################

.SUFFIXES: .f90 .F90 .c .C

%.o : %.F90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) $(INCDIRS) $(LIBDIRS) -DESMF_VERSION_MAJOR=$(ESMF_VERSION_MAJOR) $<

%.o : %.F
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) $(INCDIRS) $(LIBDIRS) -DESMF_VERSION_MAJOR=$(ESMF_VERSION_MAJOR) $<

# -----------------------------------------------------------------------------
vars:
	@echo " >>>>>>>>>>>>>>>> vars >>>>>>>>>>>>>>"
	@echo "INCDIRS $(INCDIRS)"
	@echo "LIBDIRS $(LIBDIRS)"
	@echo "ESMF_F90COMPILER $(ESMF_F90COMPILER)"
	@echo "ESMF_F90LINKER $(ESMF_F90LINKER)"
	@echo "ESMF_F90LINKOPTS $(ESMF_F90LINKOPTS)"
	@echo "ESMF_F90LINKPATHS  $(ESMF_F90LINKPATHS)"
	@echo "ESMF_F90LINKRPATHS $(ESMF_F90LINKRPATHS)"
	@echo "INCDIRS  $(INCDIRS)"
	@echo "LIBDIRS $(LIBDIRS)"
	@echo "ESMF_F90ESMFLINKLIBS $(ESMF_F90ESMFLINKLIBS)"	


# -----------------------------------------------------------------------------
esm1.x: esmApp.o esm.o wav.o adc_cap.o adc_mod.o wav_mod.o couple2swan_modif.o hwrf_mod.o hwrf_cap.o  atmesh_mod.o atmesh.o
	$(ESMF_F90LINKER) $(ESMF_F90LINKOPTS) $(ESMF_F90LINKPATHS) $(ESMF_F90LINKRPATHS) $(INCDIRS) $(LIBDIRS) -o $@ $^ $(ESMF_F90ESMFLINKLIBS) $(LIBSS) 
	echo " CP EXE FILE ..... "
	cp esm1.x  $(ADCDIR)/work/adc_esmf
    

install: $(LIBRARY)
	rm -f adcirc.mk.install
	@echo "# ESMF self-describing build dependency makefile fragment" > adcirc.mk.install
	@echo "# src location Zeus: $pwd" >> adcirc.mk.install
	@echo  >> adcirc.mk.install
	@echo "ESMF_DEP_FRONT     = adc_cap_mod" >> adcirc.mk.install
	@echo "ESMF_DEP_INCPATH   = $(INSTALLDIR)" >> adcirc.mk.install
	@echo "ESMF_DEP_CMPL_OBJS = " >> adcirc.mk.install
	@echo "ESMF_DEP_LINK_OBJS = $(ADCDIR)/work/libadc.a $(INSTALLDIR)/libadc.a" >> adcirc.mk.install
	mkdir -p $(INSTALLDIR)
	cp -f $(ADCDIR)/work/libadc.a $(INSTALLDIR)
	cp -f libadc.a adc_cap_mod.mod $(INSTALLDIR) 
	cp -f adcirc.mk.install $(INSTALLDIR)/adcirc.mk






# module dependencies:
esmApp.o:   esm.o
esm.o:      adc_cap.o
esm.o:      wav.o
esm.o:      hwrf_cap.o
esm.o:      atmesh.o
wav.o:      adc_cap.o
adc_cap.o:  couple2swan_modif.o
adc_cap.o:  adc_mod.o
wav.o:      wav_mod.o
atmesh.o:   atmesh_mod.o 
hwrf_cap.o: hwrf_mod.o


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
.PHONY: dust clean distclean info edit
dust:
	rm -f PET*.ESMF_LogFile esm1.x adc_cap.o adc_cap.mod
clean:
	rm -f esm1.x *.o *.mod *.a *.mk.install
distclean: dust clean

edit:
	gedit Makefile esmApp.F90 esm.F90 wav.F90 adc_cap.F90 adc_mod.F90 couple2swan_modif.F  hwrf_mod.F90 hwrf_cap.F90 wav_mod.F90 atmesh.F90 atmesh_mod.F90 &


