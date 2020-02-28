#!/bin/bash

## moghimis@gmail.com - 01/31/2020
## This script compiles ADCIRC model with a pre-selected sets of
## functionalities in NUOPC coupled mode in Hera computer.
# To run: ./make_nuopc esmf-impi-env.sh intel

esmf_env=$1
comp_opt=$2
compiler=intel

# environ required by ESMF/NUOPC
if [[ ! -z $esmf_env ]]; then
    echo "$0: Sourcing environ file $esmf_env"
    source $esmf_env
elif [[ ! -z $ESMFMKFILE ]]; then
    echo "$0: Setting environ for ESMFMKFILE"
else
    echo "$0: ESMFMKFILE enviorn is not set, can not continue!"
    exit 1
fi
 

if [[ ! -z $comp_opt ]]; then
    echo "$0: Setting compiler option to $comp_opt"
    export compiler=$comp_opt
else
    echo "$0: Setting compiler option by default to intel"
    export compiler=intel
fi

export NETCDFLAG=enable
export NETCDF4FLAG=enable
export NETCDF4_COMPRESSION=enable

# #make padcirc -d compiler=$compiler NETCDFLAG=enable NETCDF4FLAG=enable NETCDF4_COMPRESSION=enable
export MAKELEVEL=0
make libadc.a 

# build adcric nuopc
export ADCDIR=$(pwd)
cd ../cpl/nuopc
make -f makefile.adc_cap.nuopc nuopc





