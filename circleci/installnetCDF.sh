#!/bin/bash

installDirectory=/home/ubuntu/adcirc-cg/work/depend

if [ ! -d $installDirectory ] ; then
    mkdir $installDirectory
fi

if [ ! -d $installDirectory/netcdf_build ] ; then
    
    #...Get the netcdf-c source code
    wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.4.1.tar.gz

    #...Unzip the package
    tar -xzf netcdf-4.4.1.tar.gz

    #...Configure the code
    cd netcdf-4.4.1
    ./configure --prefix=$installDirectory/netcdf_cc_build --with-hdf5=/usr LDFLAGS="-I/usr/include -L/usr/lib"

    #...Build the code
    make V=0

    #...Install the code to the build directory
    make install

fi

#...Create links from the building directory back to
#   the directory that is cached
sudo ln -s $installDirectory/netcdf_cc_build/include/* /usr/include/.
sudo ln -s $installDirectory/netcdf_cc_build/lib/* /usr/lib/.
