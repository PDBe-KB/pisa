#!/bin/bash

srcdir=/Users/gdiazleines/programs/pisa-lite

builddir=$srcdir/build
ssmdir=$srcdir/ssm
pisalibdir=$srcdir/pisalib
mmdbdir=$srcdir/mmdb2
srsdir=$srcdir/ccp4srs
pisadir=$srcdir/pisa

cd $builddir

rm -rf $builddir/*

g++ -c -I$srcdir $mmdbdir/*.cpp
ar rvs mmdb2.a *.o
rm *.o

g++ -c -I$srcdir $ssmdir/*.cpp
ar rvs ssm.a *.o
rm *.o

g++ -c -I$srcdir $srsdir/*.cpp
ar rvs srs.a *.o
rm *.o

g++ -c -I$srcdir $pisalibdir/*.cpp
ar rvs pisalib.a *.o
rm *.o

g++ -I$srcdir $pisadir/*.cpp mmdb2.a ssm.a srs.a pisalib.a -lz -o pisa
