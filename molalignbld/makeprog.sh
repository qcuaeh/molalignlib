#!/bin/bash -e
shopt -s nullglob

if test $# -ne 2; then
   echo Error: two arguments are required
   exit 1
fi

if test ! -d "$1"; then
   echo Error: $1 does not exist
   exit 1
fi

echo Linking program...
pushd "$1" > /dev/null
"$F90" -llapack -o "$2" *.o
popd > /dev/null

echo Done
