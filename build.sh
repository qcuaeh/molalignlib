#!/bin/bash
shopt -s nullglob

# Compile source file
compile () {
   sourcefile=$sourcedir/$1
   buildfile=$build_subdir/$1
   objectfile=$build_subdir/${1%.*}.o
   objectlist+=("$objectfile")
   if ! test -e "$objectfile" || ! test -e "$buildfile" || ! diff -q "$sourcefile" "$buildfile" > /dev/null; then
      rm -f "$objectfile"
      cp -p "$sourcefile" "$buildfile"
      echo Compiling "$1"
      "$F90" "${flags[@]}" -c "$sourcefile" -o "$objectfile" -I "$build_subdir" || exit
   fi
}

repodir=$(readlink -e "$(dirname "$0")")
sourcedir=$repodir/fortran
build_dir=$repodir/__build__
bindir=$repodir/bin

set -a
. "$repodir/build.env"
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
library=false
realkind=8

options=$(getopt -o qdl:r: -- "$@") || exit
eval set -- "$options"

while true; do
   case "$1" in
   -q) quick=true; shift;;
   -d) debug=true; shift;;
   -l) library=true; libtype=$2; shift 2;;
   -r) realkind=$2; shift 2;;
   --) shift; break;;
   *) exit
   esac
done

if [[ $# -ne 0 ]]; then
   echo Error: Too many arguments
   exit
fi

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

if $library; then
   case "$libtype" in
   s)
      :
      ;;
   d)
      flags+=(-fPIC)
      ;;
   *)
      echo Invalid library type: $libtype
      exit
      ;;
   esac
fi

build_subdir=$build_dir/static/real$realkind/$comptype
so_suffix=$(grep -o '\.[a-zA-Z0-9_]*\.[a-zA-Z0-9_]*$' <<< "$(uname -r)").so

if [[ -d $build_subdir ]]; then
  if ! $quick; then
     rm -f "$build_subdir"/*
  fi
else
  mkdir -p "$build_subdir"
fi

if [[ -d $bindir ]]; then
   cd "$bindir"
   rm -f molalign
   rm -f molalignlib.a
   rm -f molalignlib${so_suffix}
   if type "$PYTHON" &> /dev/null; then
      "$PYTHON" \
<<HEREDOC
import os, sys, sysconfig
name = 'OS' if sys.version_info < (3, 4) else 'EXT_SUFFIX'
file = 'molalignlib' + sysconfig.get_config_var(name)
if os.path.exists(file): os.remove(file)
HEREDOC
   fi
else
   mkdir "$bindir"
fi

f2pylist=()
objectlist=()
cd "$build_subdir"

while IFS= read -r line; do
   eval set -- "$line"
   compile "$1"
   if [[ $2 == f2py ]]; then
      f2pylist+=("$build_subdir"/"$1")
   fi
done < <(grep -v '^#' "$sourcedir"/compilelist)

cd "$bindir"

if $library; then
   case $libtype in
   s)
      echo Linking static library...
      ar r molalignlib.a "${objectlist[@]}" &> /dev/null
      ;;
   d)
      echo Linking shared library...
      "$F90" -shared "${libpathlist[@]}" -llapack -o molalignlib${so_suffix} "${objectlist[@]}"
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
         "$F2PY" -h "$build_subdir"/molalignlib.pyf --overwrite-signature -m molalignlib "${f2pylist[@]}" --f2cmap "$f2cmap" --quiet
         "$F2PY" -c "$build_subdir"/molalignlib.pyf -I"$build_subdir" "${libpathlist[@]}" -llapack "${objectlist[@]}" --f2cmap "$f2cmap" --fcompiler=gnu95 --quiet
      fi
      ;;
   esac
else
   echo Linking program...
   "$F90" "${libpathlist[@]}" -llapack -o molalign "${objectlist[@]}"
fi
