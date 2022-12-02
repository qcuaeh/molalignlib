#!/bin/bash -e
shopt -s nullglob

clean_build() {
   if test -d build; then
      pushd build >/dev/null
      for file in *.f* *.mod *.o; do
         rm "$file"
      done
      popd >/dev/null
   fi
}

if test ! -e build.env; then
   echo Error: build.env does not exist
   exit 1
elif test ! -f build.env; then
   echo Error: build.env does exist but is not a file
   exit 1
fi

# Export environment
while IFS= read -r line; do
   var=${line%%=*}
   value=${line#*=}
   declare -- "$var"="$value"
   export -- "$var"
done < build.env

# Build program
./molalignbld/compile.sh molalignlib build
./molalignbld/compile.sh molalign build
./molalignbld/makeprog.sh build molalign
clean_build

# Run tests
./molalignbld/runtests.sh build tests/r05 -test -rec 10 -sort -fast -tol 0.17
#./molalignbld/runtests.sh build tests/r10 -test -rec 10 -sort -fast -tol 0.35
#./molalignbld/runtests.sh build tests/r20 -test -rec 10 -sort -fast -tol 0.69

# Build dynamic library
#./molalignbld/compile.sh -pic molalignlib build
#./molalignbld/makelib.sh build molalignlib.so
#clean_build

# Build extension module (python)
#./molalignbld/compile.sh -pic molalignlib build
#./molalignbld/makepyext.sh build f2py_molalignlib
#clean_build
