#!/bin/bash
shopt -s nullglob

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

quick=false
debug=false
sharedlibs=false
realkind=8

options=$(getopt -o qdlr: -- "$@") || exit
eval set -- "$options"

while true; do
   case "$1" in
   -q) quick=true; shift;;
   -d) debug=true; shift;;
   -l) sharedlibs=true; shift;;
   -r) realkind=$2; shift 2;;
   --) shift; break;;
   *) exit
   esac
done

flags=("${STAND_FLAGS[@]}")

case "$realkind" in
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
   echo Invalid precision type: $realkind
   exit
   ;;
esac

if $debug; then
   comptype=debug
   flags+=(-O0 "${DEBUG_FLAGS[@]}")
else
   comptype=standard
   flags+=(-Ofast)
fi

if $sharedlibs; then
   flags+=(-fPIC)
fi

buildir=$buildroot/static/real$realkind/$comptype
if [[ -d $buildir ]]; then
  if ! $quick; then
     rm -f "$buildir"/*
  fi
else
  mkdir -p "$buildir"
fi

if [[ -d $bindir ]]; then
   cd "$bindir"
   rm -f molalign
   rm -f molalignlib.a
   rm -f molalignlib.so
   if [[ -x $PYTHON ]]; then
      "$PYTHON" \
<<HEREDOC
import os, sys, sysconfig
name = 'OS' if sys.version_info < (3, 4) else 'EXT_SUFFIX'
file = 'molalignlib.' + sysconfig.get_config_var(name)
if os.path.exists(file): os.remove(file)
HEREDOC
   fi
else
   mkdir "$bindir"
fi

f2pylist=()
objectlist=()
cd "$buildir"

while IFS= read -r line; do
  eval set -- "$line"
  compile "$1"
  if [[ $2 == f2py ]]; then
     f2pylist+=("$buildir"/"$1")
  fi
done < <(grep -v '^#' "$sourcedir"/compilelist)

cd "$bindir"

if $sharedlibs; then
   echo Linking shared library...
   "$F90" -shared "${libpathlist[@]}" -llapack -o molalignlib.so "${objectlist[@]}"
   if [[ -n $F2PY ]]; then
      if ! type "$F2PY" &> /dev/null; then
         echo Error: F2PY is set but it was not found 
         exit
      fi
      if ! type "$PYTHON" &> /dev/null; then
         echo Error: F2PY is set but PYTHON was not found
         exit
      fi
      echo Linking shared python library...
      export PYTHONWARNINGS=ignore::Warning:setuptools.command.install
      "$F2PY" -h "$buildir"/molalignlib.pyf --overwrite-signature -m molalignlib "${f2pylist[@]}" --f2cmap "$f2cmap" --quiet
      "$F2PY" -c "$buildir"/molalignlib.pyf -I"$buildir" "${libpathlist[@]}" -llapack "${objectlist[@]}" --f2cmap "$f2cmap" \
         --fcompiler=gnu95 --quiet
   fi
else
   echo Linking program...
   "$F90" "${libpathlist[@]}" -llapack -o molalign "${objectlist[@]}"
   echo Linking static library...
   ar r molalignlib.a "${objectlist[@]}" &> /dev/null
fi
