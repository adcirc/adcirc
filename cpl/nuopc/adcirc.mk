# ESMF self-describing build dependency makefile fragment

ESMF_DEP_FRONT     = adc_cap
ESMF_DEP_INCPATH   = /scratch2/COASTAL/coastal/save/Saeed.Moghimi/models/NEMS/tests/new_nems_app/ADC-NEMS-APP-V2/ADCIRC/cpl/nuopc /scratch2/COASTAL/coastal/save/Saeed.Moghimi/models/NEMS/tests/new_nems_app/ADC-NEMS-APP-V2/ADCIRC_INSTALL 
ESMF_DEP_CMPL_OBJS = 
ESMF_DEP_LINK_OBJS =  -L/scratch2/COASTAL/coastal/save/Saeed.Moghimi/models/NEMS/tests/new_nems_app/ADC-NEMS-APP-V2/ADCIRC_INSTALL -ladc /scratch2/COASTAL/coastal/save/Saeed.Moghimi/models/NEMS/tests/new_nems_app/ADC-NEMS-APP-V2/ADCIRC_INSTALL/libadc_cap.a  -L/scratch2/COASTAL/coastal/save/Saeed.Moghimi/models/NEMS/tests/new_nems_app/ADC-NEMS-APP-V2/ADCIRC/work/  /scratch2/COASTAL/coastal/save/Saeed.Moghimi/models/NEMS/tests/new_nems_app/ADC-NEMS-APP-V2/ADCIRC/work/libadc.a  
