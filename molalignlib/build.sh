#!/bin/bash
shopt -s nullglob

# Compile source file
compile () {
   pushd "$buildir" > /dev/null
   objfile=${1%.*}.o
   objlist+=("$objfile")
   srcfile=$srcdir/$1
   if ! test -e "$objfile" || ! test -e "$1" || ! diff -q "$srcfile" "$1" > /dev/null; then
      rm -f "$objfile"
      cp -p "$srcfile" "$1"
      echo Compiling "$1"
      "$F90" "${flags[@]}" -c "$srcfile" -o "$objfile" || exit
   fi
   popd > /dev/null
}

if [[ -n $DESTDIR ]]; then
   if ! [[ -d $DESTDIR ]]; then
      echo Error: $DESTDIR does not exist or is not a directory
      exit
   fi
else
   DESTDIR=$PWD
fi

cd "$(dirname "$0")"

if [[ -f build.cfg ]]; then
   . build.cfg
else
   echo Error: build.env does not exist or is not a file
   exit
fi

if [[ -n $LAPACK_PATH ]]; then
   if [[ -d $LAPACK_PATH ]]; then
      libpathlist+=("-L$LAPACK_PATH")
   else
      echo Error: $LAPACK_PATH does not exist or is not a directory
      exit
   fi
fi

libname=molalignlib
progname=molaligncmd

quick=false
debug=false
library=false
realkind=8

options=$(getopt -o qdl:r: -- "$@") || exit
eval set -- "$options"

while true; do
   case "$1" in
   -q) quick=true; shift ;;
   -d) debug=true; shift ;;
   -l) library=true; libtype=$2; shift 2 ;;
   -r) realkind=$2; shift 2 ;;
   --) shift; break ;;
   *) exit
   esac
done

if [[ $# -ne 0 ]]; then
   echo Error: Too many arguments
   exit
fi

flags=("${BASE_FLAGS[@]}")

if $debug; then
   comptype=debug
   flags+=(-O0 "${DEBUG_FLAGS[@]}")
else
   comptype=optimized
   flags+=(-Ofast)
fi

if $library; then
   case "$libtype" in
      s) ;;
      d) flags+=(-fPIC) ;;
      py) flags+=(-fPIC) ;;
      *) echo Invalid library type: $libtype; exit
   esac
fi

case "$realkind" in
4)
   f2cmap='{"real":{"":"float"}}'
   shift
   ;;
8)
   f2cmap='{"real":{"":"double"}}'
   flags+=("${REAL8_FLAGS[@]}")
   shift
   ;;
*)
   echo Invalid precision type: $realkind
   exit
   ;;
esac

srcdir=$PWD/fortran
buildir=$DESTDIR/build/static/real$realkind/$comptype

if [[ -d $buildir ]]; then
   if ! $quick; then
      rm -f "$buildir"/*
   fi
else
   mkdir -p "$buildir"
fi

read -r -d '' python_suffix_script <<HEREDOC
import sys, sysconfig
name = 'OS' if sys.version_info < (3, 4) else 'EXT_SUFFIX'
print(sysconfig.get_config_var(name))
HEREDOC
if type "$PYTHON" &> /dev/null; then
   python_suffix=$("$PYTHON" <<< "$python_suffix_script")
fi

f2pylist=()
objlist=()

while IFS= read -r line; do
   eval set -- "$line"
   compile "$1"
   if [[ $2 == f2py ]]; then
      f2pylist+=("$1")
   fi
done < <(grep -v '^#' build.lst)

if $library; then
   case $libtype in
   s)
      echo Linking static library...
      pushd "$buildir" > /dev/null
      ar r $libname.a "${objlist[@]}" &> /dev/null
      popd > /dev/null
      mv "$buildir/$libname.a" "$DESTDIR"
      ;;
   d)
      echo Linking shared library...
      pushd "$buildir" > /dev/null
      "$F90" -shared "${libpathlist[@]}" -llapack -o $libname.so "${objlist[@]}"
      popd > /dev/null
      mv "$buildir/$libname.so" "$DESTDIR"
      ;;
   py)
      if ! type "$PYTHON" &> /dev/null; then
         echo Error: Python executable not found
         exit
      fi
      if ! type "$F2PY" &> /dev/null; then
         echo Error: F2PY executable not found 
         exit
      fi
      echo Linking shared python library...
      pushd "$buildir" > /dev/null
      echo "$f2cmap" > .f2py_f2cmap
      export PYTHONWARNINGS=ignore::Warning:setuptools.command.install
      "$F2PY" -h signature.pyf --overwrite-signature -m _$libname "${f2pylist[@]}" --quiet
      "$F2PY" -c signature.pyf "${libpathlist[@]}" -llapack "${objlist[@]}" --fcompiler=gnu95 --quiet
      popd > /dev/null
      mv "$buildir/_$libname$python_suffix" "$DESTDIR"
      ;;
   esac
else
   echo Linking program...
   pushd "$buildir" > /dev/null
   "$F90" "${libpathlist[@]}" -llapack -o $progname "${objlist[@]}"
   popd > /dev/null
   mv "$buildir/$progname" "$DESTDIR"
fi
