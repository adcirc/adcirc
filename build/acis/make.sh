#!/usr/bin/bash

## Load modules and set paths
source ./source-me-first.sh

## Make
echo 'Entering cmake'
cmake ../.. -DBUILD_ADCPREP=ON -DBUILD_ADCIRC=ON -DBUILD_ADCSWAN=ON -DBUILD_PADCIRC=ON -DBUILD_PADCSWAN=ON -DENABLE_OUTPUT_NETCDF=ON
echo 'Entering make'
make -j4 adcprep adcirc adcswan padcirc padcswan

