#!/bin/bash -e
shopt -s nullglob

if ! type "$F2PY" &>/dev/null; then
   echo Error: F2PY executable not found 
   exit 1
fi

if test $# -ne 2; then
   echo Error: two arguments are required
   exit 1
fi

if test ! -d "$1"; then
   echo Error: $1 does not exist
   exit 1
fi

echo Linking extension module...
pushd "$1" >/dev/null
"$F2PY" -h signature.pyf -m "$2" --overwrite-signature --quiet *.f2py
"$F2PY" -c signature.pyf --fcompiler=gnu95 --link-lapack --quiet *.o
popd >/dev/null

echo Done
