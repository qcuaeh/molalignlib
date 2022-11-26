#!/bin/bash -e

if test ! -e build.env; then
   echo Error: build.env does not exist
   exit 1
elif test ! -f build.env; then
   echo Error: build.env does exist but is not a file
   exit 1
fi

while IFS= read -r line; do
   var=${line%%=*}
   value=${line#*=}
   declare -- "$var"="$value"
   export -- "$var"
done < build.env

# Build program
./molalignutil/compile.sh -all molalignlib build
./molalignutil/compile.sh -all molalign build
./molalignutil/link-prog.sh build molalign

# Run tests
#./molalignutil/run-tests.sh build tests/r05 -test -rec 10 -sort -fast -tol 0.17
./molalignutil/run-tests.sh build tests/r10 -test -rec 10 -sort -fast -tol 0.35
#./molalignutil/run-tests.sh build tests/r20 -test -rec 10 -sort -fast -tol 0.69

# Build dynamic library
#./molalignutil/compile.sh -all -pic molalignlib build
#./molalignutil/compile.sh -all -pic molalign build
#./molalignutil/link-lib.sh build molalignlib.so
#./molalignutil/link-pyext.sh build f2py_molalignlib
