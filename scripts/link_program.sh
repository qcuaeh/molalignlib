#!/bin/bash -e
shopt -s nullglob

if test ! -e build.env; then
   echo Error: build.env does not exist
   exit
elif test ! -f build.env; then
   echo Error: build.env does exist but is not a file
   exit
fi

. build.env

if test $# -ne 2; then
   echo Error: two arguments are required
   exit
fi

if test ! -d "$1"; then
   echo Error: $1 does not exist
   exit
fi

echo Linking program...
pushd "$1" > /dev/null
"$F90" -llapack -o "$2" *.o
popd > /dev/null
