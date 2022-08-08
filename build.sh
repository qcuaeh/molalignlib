#!/bin/bash
shopt -s nullglob

set -a
. "$(dirname "$0")/build.env"
set +a

# Compile source file
compile () {
   sourcefile=$sourcedir/$1
   buildfile=$buildir/$1
   objectfile=$buildir/${1%.*}.o
   objectlist+=("$objectfile")
   if ! test -e "$objectfile" || ! test -e "$buildfile" || ! diff -q "$sourcefile" "$buildfile" > /dev/null; then
      rm -f "$objectfile"
      cp -p "$sourcefile" "$buildfile"
      echo Compiling "$1"
      "$F90" "${compflags[@]}" -c "$sourcefile" -o "$objectfile" -I "$buildir" -J "$buildir" || exit
   fi
}

parentdir=$(readlink -e "$(dirname "$0")")
sourcedir=$parentdir/source
buildroot=$parentdir/_build_dir
bindir=$parentdir/bin
libdir=$parentdir/lib

if [[ -n $LAPACK ]]; then
    if [[ -d $LAPACK ]]; then
        libpathlist+=("-L$LAPACK")
    else
        echo Error: Path $LAPACK does not exist or is not a directory
    fi
fi

options=$(getopt -o l:o:r:q -- "$@") || exit
eval set -- "$options"

libtype=none
realprec=double
optlevel=fast
recompile=true
compflags=()

while true; do
   case "$1" in
   -l)
      libtype=$2
      shift 2;;
   -o)
      optlevel=$2
      shift 2;;
   -r)
      realprec=$2
      shift 2;;
   -q)
      recompile=false;
      shift;;
   --)
      shift
      break;;
    *)
      exit
   esac
done

case "$realprec" in
   single) f2cmap=$sourcedir/single.f2cmap; shift;;
   double) compflags+=(-fdefault-real-8); f2cmap=$sourcedir/double.f2cmap; shift;;
   *) echo Invalid precision type: $realprec; exit; break;;
esac

case "$optlevel" in
   fast) compflags+=(-O3 -ffast-math); shift;;
   debug) compflags+=(-O0 -g -fbounds-check -fbacktrace -Wall -ffpe-trap=zero,invalid,overflow); shift;;
   *) echo Invalid optimization level: $optlevel; exit; break;;
esac

case $libtype in
none)
   buildir=$buildroot/static/$realprec/$optlevel
   if [[ -d $bindir ]]; then
      rm -f "$bindir"/molalign
   else
      mkdir "$bindir"
   fi
   ;;
static)
   buildir=$buildroot/static/$realprec/$optlevel
   if [[ -d $libdir ]]; then
      rm -f "$libdir"/molalign.a
   else
      mkdir "$libdir"
   fi
   ;;
shared)
   compflags+=(-fPIC)
   buildir=$buildroot/shared/$realprec/$optlevel
   if [[ -d $libdir ]]; then
      rm -f "$libdir"/molalign.so
   else
      mkdir "$libdir"
   fi
   ;;
python)
   compflags+=(-fPIC)
   buildir=$buildroot/shared/$realprec/$optlevel
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
   cp "$sourcedir"/molalign_wrapper.py "$libdir"
   ;;
esac
