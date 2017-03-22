#!/bin/bash

## moghimis@gmail.com



#module load intel
#module load impi
#module load netcdf
#module load hdf5
#export compiler=intel


#export NETCDFHOME=/home/Saeed.Moghimi/opt/gcc/netcdf-4.3.3.1
#export HDF5HOME=/home/Saeed.Moghimi/opt/gcc/hdf5-1.8.16-zlib

#make clean
#make clobber

export NETCDF=enable
export NETCDF4=enable
export NETCDF4_COMPRESSION=enable
export compiler=gnu

export MACHINE=x86_64
export OS=linux-gnu

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

