#!/bin/bash -e

# Build program
./molalignutil/compile.sh -all molalignlib build
./molalignutil/compile.sh -all molalign build
./molalignutil/link-prog.sh build molalign

# Build dynamic libraries
#./molalignutil/compile.sh -all -pic molalignlib build
#./molalignutil/compile.sh -all -pic molalign build
#./molalignutil/link-lib.sh build molalignlib.so
#./molalignutil/link-pylib.sh build f2py_molalignlib
