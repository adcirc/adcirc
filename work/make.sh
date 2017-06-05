#!/bin/bash

## moghimis@gmail.com

### intel compiler >>>>>>
#module load intel
#module load impi
#module load netcdf
#module load hdf5

#source /etc/profile
#module use /scratch4/NCEPDEV/nems/save/Gerhard.Theurich/Modulefiles
#module load intel impi netcdf esmf/7.0.0
export compiler=intel
###<<<<<<<<

export NETCDF=enable
export NETCDF4=enable
#export NETCDF4_COMPRESSION=enable


### GNU >>>>>
#export MACHINE=x86_64
#export compiler=gnu
#export OS=linux-gnu
### <<<<<<

#make clean
#make clobber


rm -rf libadc.a padcirc adcprep
make libadc.a 
make adcprep
#make adcirc
make padcirc
#make padcswan
#make punswan
#make hstime
#make adccmp
#make aswip
#make adcpost
#make adcprep_be
#make adcswan
#make p15
#make owi22
#make build13
#make buildstwave23
#make hot2asc
#make inflate

