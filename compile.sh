#!/bin/bash

# Modify srcdir below to point to complete path of the PISA-LITE directory
# For example: srcdir=/Users/foo/bar/pisa-lite

# srcdir=/Users/gdiazleines/programs/pisa-lite
# srcdir=/media/jellaway/FlashData/EMBL-EBI/pisa_work/pisa-lite
srcdir=$SRCDIR

mkdir build

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
g++ -I$srcdir $pisadir/*.cpp *.o -lz -o pisa
