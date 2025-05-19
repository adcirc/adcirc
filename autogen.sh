#!/usr/bin/env bash

system_type=$(uname -s)
if [ x"$system_type" == "xDarwin" ]; then
  LIBTOOL=glibtoolize
else
  LIBTOOL=libtoolize
fi  

$LIBTOOL --copy
aclocal
automake --add-missing --copy
autoconf
