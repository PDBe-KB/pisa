#!/bin/bash
set -e

# Modify srcdir below to point to complete path of the PISA directory                                                                              
# For example: srcdir=/Users/foo/bar/pisa
# Ensure SRCDIR is set
if [ -z "$SRCDIR" ]; then
    echo "Error: SRCDIR is not set"
    exit 1
fi

srcdir="$SRCDIR"

builddir="$srcdir/build"
ssmdir="$srcdir/ssm"
pisalibdir="$srcdir/pisalib"
mmdbdir="$srcdir/mmdb2"
srsdir="$srcdir/ccp4srs"
pisadir="$srcdir/pisa"

mkdir -p "$builddir"
cd "$builddir"

# compile each component
g++ -c -I"$srcdir" "$mmdbdir"/*.cpp -Wno-deprecated-declarations
g++ -c -I"$srcdir" "$ssmdir"/*.cpp -Wno-deprecated-declarations
g++ -c -I"$srcdir" "$srsdir"/*.cpp -Wno-deprecated-declarations
g++ -c -I"$srcdir" "$pisalibdir"/*.cpp -Wno-deprecated-declarations
g++ -c -I"$srcdir" "$pisadir"/*.cpp -Wno-deprecated-declarations

# link all objects
g++ *.o -lz -o pisa
