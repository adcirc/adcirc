#!/usr/bin/bash

## Load modules
module load intelc/18.0.0
module load intelfort/18.0.0
MODULEPATH=$MODULEPATH:/projects/acis/modules/modulefiles
module load hdf5/1.8.12-acis
module load netcdf/4.1.2-acis
module load mvapich2/2.0-acis
module load zlib/1.2.11_intel-18.0.0
module list

## Set paths
export NETCDFHOME=`nc-config --prefix`
echo 'NETCDFHOME is set to ':
echo "  $NETCDFHOME"

