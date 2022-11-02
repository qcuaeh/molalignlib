#!/bin/bash
shopt -s nullglob

# Compile source file
compile () {
   cd "$buildir"
   sourcefile=$sourcedir/$1
   buildfile=$buildir/$1
   objectfile=$buildir/${1%.*}.o
   objectlist+=("$objectfile")
   if ! test -e "$objectfile" || ! test -e "$buildfile" || ! diff -q "$sourcefile" "$buildfile" > /dev/null; then
      rm -f "$objectfile"
      cp -p "$sourcefile" "$buildfile"
      echo Compiling "$1"
      "$F90" "${flags[@]}" -c "$sourcefile" -o "$objectfile" -I "$buildir" || exit
   fi
}

parentdir=$(readlink -e "$(dirname "$0")")
sourcedir=$parentdir/source
buildroot=$parentdir/__build__
bindir=$parentdir/bin
libdir=$parentdir/lib

set -a
. "$parentdir/build.env"
set +a

if [[ -n $LAPACK_PATH ]]; then
    if [[ -d $LAPACK_PATH ]]; then
        libpathlist+=("-L$LAPACK_PATH")
    else
        echo Error: Path $LAPACK_PATH does not exist or is not a directory
    fi
fi

realprec=8
libtype=none
comptype=fast
recompile=true

options=$(getopt -o l:r:dq -- "$@") || exit
eval set -- "$options"

while true; do
   case "$1" in
   -l) libtype=$2; shift 2;;
   -r) realprec=$2; shift 2;;
   -d) comptype=debug; shift;;
   -q) recompile=false; shift;;
   --) shift; break;;
   *) exit
   esac
done

flags=("${STAND_FLAGS[@]}")

case "$realprec" in
4)
   f2cmap=$sourcedir/f2cmap_single
   shift
   ;;
8)
   f2cmap=$sourcedir/f2cmap_double
   flags+=("${REAL8_FLAGS[@]}")
   shift
   ;;
*)
   echo Invalid precision type: $realprec
   exit
   ;;
esac

case "$comptype" in
   fast) flags+=(-Ofast); shift;;
   debug) flags+=(-O0 "${DEBUG_FLAGS[@]}"); shift;;
   *) echo Invalid optimization level: $comptype; exit;;
esac

case $libtype in
none)
   buildir=$buildroot/static/real$realprec/$comptype
   if [[ -d $bindir ]]; then
      rm -f "$bindir"/molalign
   else
      mkdir "$bindir"
   fi
   ;;
static)
   buildir=$buildroot/static/real$realprec/$comptype
   if [[ -d $libdir ]]; then
      rm -f "$libdir"/molalign.a
   else
      mkdir "$libdir"
   fi
   ;;
shared)
   flags+=(-fPIC)
   buildir=$buildroot/dynamic/real$realprec/$comptype
   if [[ -d $libdir ]]; then
      rm -f "$libdir"/molalign.so
   else
      mkdir "$libdir"
   fi
   ;;
python)
   flags+=(-fPIC)
   buildir=$buildroot/dynamic/real$realprec/$comptype
   if [[ -d $libdir ]]; then
      "$PYTHON" \
<<HEREDOC
import os, sys, sysconfig
name = 'OS' if sys.version_info < (3, 4) else 'EXT_SUFFIX'
file = 'molalign.' + sysconfig.get_config_var(name)
if os.path.exists(file): os.remove(file)
HEREDOC
   else
      mkdir "$libdir"
   fi
   ;;
*)
   echo Invalid library type: $libtype
   exit
esac

if [[ -d $buildir ]]; then
  if $recompile; then
     rm -f "$buildir"/*
  fi
else
  mkdir -p "$buildir"
fi

objectlist=()
f2pylist=()

while IFS= read -r line; do
  eval set -- "$line"
  compile "$1"
  if [[ $2 == f2py ]]; then
     f2pylist+=("$buildir"/"$1")
  fi
done < <(grep -v '^#' "$sourcedir"/compilelist)

case $libtype in
none)
   cd "$bindir"
   echo Linking program...
   "$F90" "${libpathlist[@]}" -llapack -o molalign "${objectlist[@]}"
   ;;
static)
   cd "$libdir"
   echo Linking static library...
   ar r molalign.a "${objectlist[@]}"
   ;;
shared)
   cd "$libdir"
   echo Linking shared library...
   "$F90" -shared "${libpathlist[@]}" -llapack -o molalign.so "${objectlist[@]}"
   ;;
python)
   cd "$libdir"
   echo Linking python library...
   export PYTHONWARNINGS=ignore::Warning:setuptools.command.install
   "$F2PY" -h "$buildir"/molalign.pyf --overwrite-signature -m molalign "${f2pylist[@]}" --f2cmap "$f2cmap" --quiet
   "$F2PY" -c "$buildir"/molalign.pyf -I"$buildir" "${libpathlist[@]}" -llapack "${objectlist[@]}" --f2cmap "$f2cmap" \
      --fcompiler=gnu95 --quiet
   cp "$sourcedir"/wrapper.py "$libdir"
   ;;
esac
