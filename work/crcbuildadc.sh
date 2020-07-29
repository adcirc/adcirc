#!/bin/bash

###########################################
# exaple usages: 
#    make padcirc
#    make adcprep 
#    make adcprep padcirc
###########################################

# NETCDF locations
export NETCDF=/opt/crc/n/netcdf/4.7.0/intel/18.0/      
export NETCDFHOME=/opt/crc/n/netcdf/4.7.0/intel/18.0/
# export ESMFMKFILE=/afs/crc.nd.edu/user/d/dwirasae/AlaskaProject/esmf_8/DEFAULTINSTALLDIR/lib/libO/Linux.intel.64.mvapich2.default/esmf.mk 

# compiler vendor 
cpler=intel 

for var in  "$@"
do
 echo "Building "$var" ......" 
 make compiler=$cpler NETCDF=enable NETCDF4=enable $var 
done



