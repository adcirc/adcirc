# ESMF self-describing build dependency makefile fragment

ESMF_DEP_FRONT     = adc_cap
ESMF_DEP_INCPATH   = /scratch2/COASTAL/coastal/save/shared/repositories/ADC-WW3-NWM-NEMS/ADCIRC/thirdparty/nuopc /scratch2/COASTAL/coastal/save/shared/repositories/ADC-WW3-NWM-NEMS/ADCIRC_INSTALL 
ESMF_DEP_CMPL_OBJS = 
ESMF_DEP_LINK_OBJS =  -L/scratch2/COASTAL/coastal/save/shared/repositories/ADC-WW3-NWM-NEMS/ADCIRC_INSTALL -ladc /scratch2/COASTAL/coastal/save/shared/repositories/ADC-WW3-NWM-NEMS/ADCIRC_INSTALL/libadc_cap.a  -L/scratch2/COASTAL/coastal/save/shared/repositories/ADC-WW3-NWM-NEMS/ADCIRC/work/  /scratch2/COASTAL/coastal/save/shared/repositories/ADC-WW3-NWM-NEMS/ADCIRC/work/libadc.a  
