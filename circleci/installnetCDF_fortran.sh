#!/bin/bash

installDirectory=/home/ubuntu/adcirc-cg/work/depend

if [ ! -d $installDirectory ] ; then
    mkdir $installDirectory
fi

if [ ! -d $installDirectory/netcdf_fortran_build ] ; then
    
    #...Get the netcdf-c source code
    wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-fortran-4.4.4.tar.gz

    #...Unzip the package
    tar -xzf netcdf-fortran-4.4.4.tar.gz

    #...Configure the code
    cd netcdf-fortran-4.4.4
    ./configure --prefix=$installDirectory/netcdf_fortran_build LDFLAGS="-I/usr/include -L/usr/lib"

    #...Build the code
    make V=0

    #...Install the code to the build directory
    make install

fi

#...Create links from the building directory back to
#   the directory that is cached
sudo ln -s $installDirectory/netcdf_fortran_build/include/* /usr/include/.
sudo ln -s $installDirectory/netcdf_fortran_build/lib/* /usr/lib/.
