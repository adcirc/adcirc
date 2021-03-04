# ESMF self-describing build dependency makefile fragment

ESMF_DEP_FRONT     = adc_cap
ESMF_DEP_INCPATH   = /work/07531/zrb/stampede2/nems/builds/ADC-WW3-NWM-NEMS/ADCIRC/thirdparty/nuopc /work/07531/zrb/stampede2/nems/builds/ADC-WW3-NWM-NEMS/ADCIRC_INSTALL 
ESMF_DEP_CMPL_OBJS = 
ESMF_DEP_LINK_OBJS =  -L/work/07531/zrb/stampede2/nems/builds/ADC-WW3-NWM-NEMS/ADCIRC_INSTALL -ladc /work/07531/zrb/stampede2/nems/builds/ADC-WW3-NWM-NEMS/ADCIRC_INSTALL/libadc_cap.a  -L/work/07531/zrb/stampede2/nems/builds/ADC-WW3-NWM-NEMS/ADCIRC/work/  /work/07531/zrb/stampede2/nems/builds/ADC-WW3-NWM-NEMS/ADCIRC/work/libadc.a  
