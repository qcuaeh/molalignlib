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

if ! type "$F2PY" &>/dev/null; then
   echo Error: F2PY executable not found 
   exit
fi

if test $# -ne 2; then
   echo Error: two arguments are required
   exit
fi

if test ! -d "$1"; then
   echo Error: $1 does not exist
   exit
fi

PYTHONWARNINGS=ignore::Warning:setuptools.command.install
export PYTHONWARNINGS

pushd "$1" >/dev/null
echo Linking dynamic python library...
"$F2PY" -h signature.pyf -m "$2" --overwrite-signature --quiet *.f2py
"$F2PY" -c signature.pyf --fcompiler=gnu95 -llapack --quiet *.o
popd >/dev/null
