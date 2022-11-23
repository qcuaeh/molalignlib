#!/bin/bash -e

# Build program
./scripts/compile.sh -all molalignlib build
./scripts/compile.sh -all molalign build
./scripts/link_program.sh build molalign

# Build dynamic libraries
#./scripts/compile.sh -all -pic molalignlib build
#./scripts/compile.sh -all -pic molalign build
#./scripts/link_library.sh build molalignlib.so
#./scripts/link_pythonlib.sh build f2py_molalignlib
