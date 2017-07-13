#!/bin/bash

## moghimis@gmail.com
## ADCIRC  for Theia

export compiler=intel
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


rm -rf libadc.a #padcirc adcprep
make libadc.a 
make adcprep
#make adcirc
#make padcirc
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

