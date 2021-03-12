#!/bin/bash

# Panagiotis Velissariou <panagiotis.velissariou@noaa.gov> - 12/05/2020
# Original by: moghimis@gmail.com - 01/31/2020
# This script compiles ADCIRC model with a pre-selected sets of
## functionalities in NUOPC coupled mode.
# Run it as: ./make_nuopc.sh some-env-module-file some-compiler-name   or,
#            ./make_nuopc.sh some-compiler-name                        or,
#            ./make_nuopc.sh
# Note: some-env-module-file will be sourced from this bash script
#       so it should be bash compatible file.

# Get the directory where the script is located
if [[ $(uname -s) == Darwin ]]; then
  readonly scrNAME="$(basename "$(grealpath -s "${BASH_SOURCE[0]}")")"
  readonly scrDIR="$(cd "$(dirname "$(grealpath -s "${BASH_SOURCE[0]}")")" && pwd -P)"
else
  readonly scrNAME="$(basename "$(realpath -s "${BASH_SOURCE[0]}")")"
  readonly scrDIR="$(cd "$(dirname "$(realpath -s "${BASH_SOURCE[0]}")")" && pwd -P)"
fi

###====================
### Get the script arguments: env. file and compiler
###====================
num_opt=${#}
if [ ${num_opt} -ge 2 ]; then
  esmf_env="${1}"
  comp_opt="$(echo "${2}" | tr '[:upper:]' '[:lower:]')"
elif [ ${num_opt} -eq 1 ]; then
  comp_opt="$(echo "${1}" | tr '[:upper:]' '[:lower:]')"
else
  if [ -n ${NEMS_COMPILER:+1} ]; then
    comp_opt="$(echo "${NEMS_COMPILER}" | tr '[:upper:]' '[:lower:]')"
  fi
fi

# Environment required by ESMF/NUOPC (OPTIONAL)
if [ -n "${esmf_env:+1}" ]; then
  if [ -f "${esmf_env}" ]; then
    echo "${scrNAME} :: Sourcing the environment file: \"${esmf_env}\""
    source ${esmf_env}
  else
    echo "${scrNAME} :: Using the the environment variable ESMFMKFILE"
  fi
fi

# Compiler to use (REQUIRED)
case "${comp_opt}" in
  gnu   ) ;; #Do nothing
  intel ) ;; #Do nothing
  pgi   ) ;; #Do nothing
  *     ) echo "${scrNAME} :: Compiler \"${comp_opt}\" is not supported"
          echo "   Exiting ..."
          exit 1
          ;;
esac
###====================

###====================
### Set/Export the important environment variables
###====================
if [[ -n ${NETCDFHOME:+1} ]]; then
  export NETCDFLAG=enable
  export NETCDF4FLAG=enable
  export NETCDF4_COMPRESSION=enable
else
  echo "${scrNAME} :: NETCDFHOME is not defined. Define this environment variable before running this script."
  echo "   Exiting ..."
  exit 1
fi

export compiler=${compiler}
export ADCDIR=$(dirname $(dirname ${scrDIR}))
###====================


# move to the `work/` directory to build main ADCIRC executables
pushd ${ADCDIR}/work >/dev/null 2>&1
  export MAKELEVEL=0
  make NETCDF=enable NETCDF4=enable libadc.a
  make NETCDF=enable NETCDF4=enable adcprep
  #make padcirc -d compiler=$compiler NETCDFLAG=enable NETCDF4FLAG=enable NETCDF4_COMPRESSION=enable
popd >/dev/null 2>&1

# build ADCIRC NUOPC
make -f makefile.adc_cap.nuopc nuopc

# build ESTOFS `tidefac` executable
pushd ${ADCDIR}/work >/dev/null 2>&1
  if [ ! -d util ]; then
    mkdir -p util
  else
    [ -f util/tidefac ] && rm -f util/tidefac
  fi
popd >/dev/null 2>&1
pushd ${ADCDIR}/util/estofs_tide_fac >/dev/null 2>&1
  make
  cp -f estofs_tide_fac ${ADCDIR}/work/util/tidefac
popd >/dev/null 2>&1
