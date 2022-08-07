#!/bin/bash
shopt -s nullglob

set -a
. "$(dirname "$0")/build.env"
set +a

# Compile source file
compile () {
   sourcefile=$sourcedir/$1
   buildfile=$builddir/$1
   objectfile=$builddir/${1%.*}.o
   objectlist+=("$objectfile")
   if ! test -e "$objectfile" || ! test -e "$buildfile" || ! diff -q "$sourcefile" "$buildfile" > /dev/null; then
      rm -f "$objectfile"
      cp -p "$sourcefile" "$buildfile"
      echo Compiling "$1"
      "$F90" "${compflags[@]}" -c "$sourcefile" -o "$objectfile" -I "$builddir" -J "$builddir" || exit
   fi
}

parentdir=$(readlink -e "$(dirname "$0")")
sourcedir=$parentdir/source
buildroot=$parentdir/_build_dir
bindir=$parentdir/bin

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
   builddir=$buildroot/static/$realprec/$optlevel
   if [[ -d $bindir ]]; then
      rm -f "$bindir"/molalign
   else
      mkdir "$bindir"
   fi
   ;;
static)
   builddir=$buildroot/static/$realprec/$optlevel
   if [[ -d $bindir ]]; then
      rm -f "$bindir"/molalign.a
   else
      mkdir "$bindir"
   fi
   ;;
shared)
   compflags+=(-fPIC)
   builddir=$buildroot/shared/$realprec/$optlevel
   if [[ -d $bindir ]]; then
      rm -f "$bindir"/molalign.so
   else
      mkdir "$bindir"
   fi
   ;;
python)
   compflags+=(-fPIC)
   builddir=$buildroot/shared/$realprec/$optlevel
   if [[ -d $bindir ]]; then
      "$PYTHON" \
<<HEREDOC
import os, sys, sysconfig
name = 'OS' if sys.version_info < (3, 3, 1) else 'EXT_SUFFIX'
file = 'molalign.' + sysconfig.get_config_var(name)
if os.path.exists(file): os.remove(file)
HEREDOC
   else
      mkdir "$bindir"
   fi
   ;;
*)
   echo Invalid library type: $libtype
   exit
esac

if [[ -d $builddir ]]; then
  if $recompile; then
     rm -f "$builddir"/*
  fi
else
  mkdir -p "$builddir"
fi

objectlist=()
exportlist=()

while IFS= read -r line; do
  eval set -- "$line"
  compile "$1"
  if [[ -n $2 ]]; then
     exportlist+=("$builddir"/"$1")
  fi
done < <(grep -v '^#' "$sourcedir"/compilelist)

cd "$bindir"

case $libtype in
none)
   echo Linking program...
   "$F90" "${libpathlist[@]}" -llapack -o molalign "${objectlist[@]}"
   ;;
static)
   echo Linking static library...
   ar r molalign.a "${objectlist[@]}"
   ;;
shared)
   echo Linking shared library...
   "$F90" -shared "${libpathlist[@]}" -llapack -o molalign.so "${objectlist[@]}"
   ;;
python)
   echo Linking python library...
   export PYTHONWARNINGS=ignore::Warning:setuptools.command.install
   "$F2PY" -h "$builddir"/molalign.pyf --overwrite-signature -m molalign "${exportlist[@]}" --f2cmap "$f2cmap" --quiet
   "$F2PY" -c "$builddir"/molalign.pyf -I"$builddir" "${libpathlist[@]}" -llapack "${objectlist[@]}" --f2cmap "$f2cmap" \
      --fcompiler=gnu95 --quiet
   ;;
esac
