#!/bin/bash

installDirectory=/home/ubuntu/adcirc-cg/work/depend

if [ ! -d $installDirectory/xdmf_build ] ; then

    if [ ! -s $installDirectory ] ; then
        mkdir $installDirectory
    fi

    cd $installDirectory
    
    #...Get the xdmf source code
    git clone https://gitlab.kitware.com/xdmf/xdmf.git

    #...Configure the code
    mkdir xdmf/build
    cd xdmf/build
    cmake -DCMAKE_BUILD_TYPE=Release -DXDMF_BUILD_UTILS=ON -DCMAKE_INSTALL_PREFIX=../../xdmf_build -DXDMF_BUILD_FORTRAN=ON -DBUILD_SHARED_LIBS=ON -DCMAKE_SHARED_LINKER_FLAGS=-fPIC -DCMAKE_CXX_FLAGS=-fPIC  ../

    #...Build the code using 4 processors
    make

    #...Install the code to the build directory
    make install

fi

#...Create links from the building directory back to
#   the directory that is cached
sudo ln -sf $installDirectory/xdmf_build/include/* /usr/include/.
sudo ln -sf $installDirectory/xdmf_build/lib/lib/*.so /usr/lib/.
sudo ln -sf $installDirectory/xdmf_build/lib/x86_64-linux-gnu/*.so /usr/lib/.
