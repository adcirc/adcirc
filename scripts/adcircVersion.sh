#!/bin/bash

#............................#
#     adcircVersion          #
#                            #
# Script for automatically   #
# creating version.F         #
#                            #
#............................#

#...Check if called with correct arguments
if [ $# != 1 ] ; then
    echo "ERROR: Run with adcirc root directory"
    exit 1
fi

dir=$1

#...Get the version from version.F as our default
if [ -s $dir/version_default.F ] ; then
    defaultVersion=$(cat $dir/version_default.F | grep ADC_VERSION | cut -d= -f2)
    defaultVersion=$(echo $defaultVersion) #...Strip leading chars
    defaultVersion="${defaultVersion%\"}"  #...Strip leading "
    defaultVersion="${defaultVersion#\"}"  #...Strip trailing "
    foundDefaultVersion=1
fi

#...Get the version from git. If we are not
#   within a repository (signaled by a git error)
#   quit and only write what is already 
#   contained in version.F
version=$(git describe --always --tags 2>/dev/null)
if [ "x$?" != "x0" ] ; then
    if [ ! -s $dir/version.F ] ; then
        echo $defaultVersion
        cat $dir/version_default.F > $dir/version.F
        echo $defaultVersion
    else
        foundVersion=$(cat $dir/version.F | grep ADC_VERSION | cut -d= -f2)
        foundVersion=$(echo $foundVersion) #...Strip leading chars
        foundVersion="${foundVersion%\"}"  #...Strip leading "
        foundVersion="${foundVersion#\"}"  #...Strip trailing "
        echo $foundVersion 
    fi
    exit 0
fi

#...Get the full hash for the current state of the
#   repository
adcirc_hash=$(git rev-parse HEAD)

#...If version_autogen.F exists, check to make sure
#   we don't trigger a rebuild by rewriting it
if [ -s $dir/version.F ] ; then
    versionFile=$(cat $dir/version.F | grep ADC_VERSION)
    versionGit="character(80), parameter :: ADC_VERSION = \""$version\"
    if [ "x$versionLine" == "x$versionGit" ] ; then
        echo $version
        exit 0
    fi
fi

#...Write an updated version.F with the descriptive commit id and echo
#   the version for the makefile
echo "      module version" > $dir/version.F
echo "      character(80), parameter :: ADC_VERSION = \""$version\" >> $dir/version.F
echo "      character(80), parameter :: ADC_HASH    = \""$adcirc_hash\" >> $dir/version.F
echo "      integer, save       :: FileFmtMajor = 1" >> $dir/version.F
echo "      integer, save       :: FileFmtMinor = 2" >> $dir/version.F
echo "      integer, save       :: FileFmtRev   = 0" >> $dir/version.F
echo "      end module version" >> $dir/version.F
echo $version
exit 0
