#!/usr/bin/env bash
libtoolize --copy
aclocal
automake --add-missing --copy
autoconf
