#!/bin/bash

# Panagiotis Velissariou <panagiotis.velissariou@noaa.gov> - 07/20/2021
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
  if [ -n "${NEMS_COMPILER:+1}" ]; then
    comp_opt="$(echo "${NEMS_COMPILER}" | tr '[:upper:]' '[:lower:]')"
  fi
fi

# Environment required by ESMF/NUOPC (OPTIONAL)
if [ -n "${esmf_env:+1}" ]; then
  if [ -f "${esmf_env}" ]; then
    echo "${scrNAME} :: Sourcing the environment file: \"${esmf_env}\""
    source ${esmf_env}
  else
    echo "${scrNAME} :: The environment file: \"${esmf_env}\" does not exist"
    echo "${scrNAME} :: Trying the environment variable ESMFMKFILE ..."
    if [ -n "${ESMFMKFILE:+1}" ]; then
      echo "${scrNAME} :: Using the environment variable ESMFMKFILE"
    else
      echo "${scrNAME} :: WARNING: The environment variable ESMFMKFILE is not defined"
    fi
  fi
else
  if [ -n "${ESMFMKFILE:+1}" ]; then
    echo "${scrNAME} :: Using the environment variable ESMFMKFILE"
  else
    echo "${scrNAME} :: WARNING: The environment variable ESMFMKFILE is not defined"
  fi
fi

# Compiler to use (REQUIRED)
case "${comp_opt}" in
  gnu   ) ;; #Do nothing
  intel ) ;; #Do nothing
  pgi   ) ;; #Do nothing
  *     ) echo "${scrNAME} :: Compiler \"${comp_opt:-UNDEF}\" is not supported"
          echo "   Exiting ..."
          exit 1
          ;;
esac

adc_exe="$(echo "${BUILD_EXECS:-UNDEF}" | tr '[:upper:]' '[:lower:]')"
###====================

###====================
### Set/Export the important environment variables
###====================
if [[ -n ${NETCDFHOME:+1} ]]; then
  export NETCDF=enable
  export NETCDF4=enable
  export NETCDF4_COMPRESSION=enable
else
  echo "${scrNAME} :: NETCDFHOME is not defined. Define this environment variable before running this script."
  echo "   Exiting ..."
  exit 1
fi

export ADCDIR=$(dirname $(dirname ${scrDIR}))
###====================

# move to the `work/` directory to build main ADCIRC executables
pushd ${ADCDIR}/work >/dev/null 2>&1
  export MAKELEVEL=0

  # Mandatory components
  make compiler=${comp_opt} libadc.a
  make compiler=${comp_opt} adcprep

  # Optional components
  for iexe in ${adc_exe}
  do
    case "${iexe}" in
      adcirc)
        make compiler=${comp_opt} adcirc
        ;;
      padcirc)
        make compiler=${comp_opt} padcirc
        ;;
      aswip)
        make compiler=${comp_opt} aswip
        ;;
      *) ;; #Do nothing
    esac
  done

popd >/dev/null 2>&1

# build ADCIRC NUOPC
make nuopc

# build `tide_fac` executable from ESTOFS
pushd ${ADCDIR}/work >/dev/null 2>&1
  if [ ! -d util ]; then
    mkdir -p util
  else
    [ -f util/tide_fac ] && rm -f util/tide_fac
  fi
popd >/dev/null 2>&1
pushd ${ADCDIR}/util/estofs_tide_fac >/dev/null 2>&1
  make
  cp -f estofs_tide_fac ${ADCDIR}/work/util/tide_fac
popd >/dev/null 2>&1
