#!/usr/bin/env bash
libtoolize
aclocal
automake --add-missing
autoconf
