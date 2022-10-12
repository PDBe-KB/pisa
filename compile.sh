#!/bin/bash

srcdir=$PWD
builddir=$srcdir/build
ssmdir=$srcdir/ssm
pisalibdir=$srcdir/pisalib
mmdbdir=$srcdir/mmdb2
srsdir=$srcdir/ccp4srs
pisadir=$srcdir/pisa

cd $builddir


g++ -c -I$srcdir $mmdbdir/*.cpp

g++ -c -I$srcdir $ssmdir/*.cpp

g++ -c -I$srcdir $srsdir/*.cpp

g++ -c -I$srcdir $pisalibdir/*.cpp

# Weslley changes
g++ -I$srcdir $pisadir/*.cpp *.o -lz -o pisa



