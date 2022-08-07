#!/bin/bash
shopt -s nullglob
parentdir=$(cd -- "$(dirname "$0")" && pwd)
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

sourcedir=$parentdir/source
buildroot=$parentdir/_build_dir
testdir=$parentdir/tests

if [[ -n $LAPACK ]]; then
    if [[ -d $LAPACK ]]; then
        libpathlist+=("-L$LAPACK")
    else
        echo Error: Path $LAPACK does not exist or is not a directory
    fi
fi

options=$(getopt -a -o '' -l program,library,slow,fast,debug,single,double,quick -- "$@") || exit
eval set -- "$options"

buildtype=program
realprec=double
optlevel=optimized
recompile=true

while true; do
   case "$1" in
   --program) buildtype=program; shift;;
   --library) buildtype=library; shift;;
   --optimized) optlevel=optimized; shift;;
   --debug) optlevel=debug; shift;;
   --single) realprec=single; shift;;
   --double) realprec=double; shift;;
   --quick) recompile=false; shift;;
   --) shift; break;;
   esac
done

builddir=$buildroot/$buildtype/$realprec/$optlevel

case $buildtype in
   program) :;;
   library) compflags+=(-fPIC);;
   *) echo Invalid build type: $buildtype; exit;;
esac

case "$optlevel" in
   optimized) compflags+=(-O3 -ffast-math); shift;;
   debug) compflags+=(-O0 -g -fbounds-check -fbacktrace -Wall -ffpe-trap=zero,invalid,overflow); shift;;
   *) echo Invalid build type: $optlevel; exit; break;;
esac

case "$realprec" in
   single) f2cmap=$sourcedir/single.f2cmap shift;;
   double) compflags+=(-fdefault-real-8); f2cmap=$sourcedir/double.f2cmap; shift;;
   *) echo Invalid precision type: $realprec; exit; break;;
esac

case $buildtype in
program)
   if [[ -d $testdir ]]; then
      rm -f "$testdir"/molalign
   else
      mkdir "$testdir"
   fi
   ;;
library)
   if [[ -d $testdir ]]; then
      rm -f "$testdir"/molalign.*.so
   else
      mkdir "$testdir"
   fi
   ;;
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

case $buildtype in
program)
   echo Linking program...
   "$F90" "${libpathlist[@]}" -llapack -o "$testdir"/molalign "${objectlist[@]}"
   ;;
library)
   cd "$testdir"
   echo Linking libraries...
   export PYTHONWARNINGS=ignore::Warning:setuptools.command.install
   "$F90" -shared "${libpathlist[@]}" -llapack -o "$testdir"/molalign.$(uname -i)-$(uname -s).so "${objectlist[@]}"
   "$F2PY" -h "$builddir"/molalign.pyf --overwrite-signature -m molalign "${exportlist[@]}" --f2cmap "$f2cmap" --quiet
   "$F2PY" -c "$builddir"/molalign.pyf -I"$builddir" "${libpathlist[@]}" -llapack "${objectlist[@]}" --f2cmap "$f2cmap" \
       --fcompiler=gnu95 --quiet
esac
