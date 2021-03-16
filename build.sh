#!/bin/bash
#This is a build script that helps with the initial build as described in README file
libtoolize --force
aclocal
autoheader
automake -ac
autoconf
./configure
make clean
make
sudo make install
echo "Warning!"
echo "Python scripts for running the trisurf-ng have been moved to separate package"

