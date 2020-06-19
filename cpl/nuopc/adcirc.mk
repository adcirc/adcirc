# ESMF self-describing build dependency makefile fragment

ESMF_DEP_FRONT     = adc_cap
ESMF_DEP_INCPATH   = /afs/crc.nd.edu/user/d/dwirasae/AlaskaProject/SaeeedNEM/ADC-WW3-NWM-NEMS/ADCIRC-CG-v55/cpl/nuopc /afs/crc.nd.edu/user/d/dwirasae/AlaskaProject/SaeeedNEM/ADC-WW3-NWM-NEMS/ADCIRC_INSTALL 
ESMF_DEP_CMPL_OBJS = 
ESMF_DEP_LINK_OBJS =  -L/afs/crc.nd.edu/user/d/dwirasae/AlaskaProject/SaeeedNEM/ADC-WW3-NWM-NEMS/ADCIRC_INSTALL -ladc /afs/crc.nd.edu/user/d/dwirasae/AlaskaProject/SaeeedNEM/ADC-WW3-NWM-NEMS/ADCIRC_INSTALL/libadc_cap.a  -L/afs/crc.nd.edu/user/d/dwirasae/AlaskaProject/SaeeedNEM/ADC-WW3-NWM-NEMS/ADCIRC/work/  /afs/crc.nd.edu/user/d/dwirasae/AlaskaProject/SaeeedNEM/ADC-WW3-NWM-NEMS/ADCIRC/work/libadc.a  
