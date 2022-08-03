# Compile source file
compile () {
   sourcefile=$SRCDIR/$1
   buildfile=$BUILDIR/$1
   objectfile=$BUILDIR/${1%.*}.o
   objectlist+=("$objectfile")
   if ! test -e "$objectfile" || ! test -e "$buildfile" || ! diff -q "$sourcefile" "$buildfile" > /dev/null; then
      rm -f "$objectfile"
      cp -p "$sourcefile" "$buildfile"
      echo Compiling "$1"
      $FORTRAN "${compflags[@]}" -c "$sourcefile" -o "$objectfile" -I "$BUILDIR" -J "$BUILDIR" || exit
   fi
}

NAME=ralign
HERE=$(cd -- "$(dirname "$0")" && pwd)
SRCDIR=$HERE/source
BINDIR=$HERE/bin
BUILDROOT=$HERE/_build_dir

if [[ -n $LAPACK_LIBRARY_PATH ]]; then
    if [[ -d $LAPACK_LIBRARY_PATH ]]; then
        linkflags+=("-L$LAPACK_LIBRARY_PATH")
    else
        echo Error: Path $LAPACK_LIBRARY_PATH does not exist or is not a directory
    fi
fi

shopt -s nullglob

options=$(getopt -a -o '' -l program,library,slow,fast,debug,single,double,recompile -- "$@") || exit
eval set -- "$options"

buildtype=program
realprec=double
optlevel=fast
recompile=false

while true; do
   case "$1" in
   --program) buildtype=program; shift;;
   --library) buildtype=library; shift;;
   --slow) optlevel=slow; shift;;
   --fast) optlevel=fast; shift;;
   --debug) optlevel=debug; shift;;
   --single) realprec=single; shift;;
   --double) realprec=double; shift;;
   --recompile) recompile=true; shift;;
   --) shift; break;;
   esac
done

BUILDIR=$BUILDROOT/$buildtype/$realprec/$optlevel

case $buildtype in
   program) :;;
   library) compflags+=(-fPIC);;
   *) echo Invalid build type: $buildtype; exit;;
esac

case "$optlevel" in
   slow) compflags+=(-O0); shift;;
   fast) compflags+=(-O3 -ffast-math); shift;;
   debug) compflags+=(-O0 -g -fbounds-check -fbacktrace -Wall -ffpe-trap=zero,invalid,overflow); shift;;
   *) echo Invalid build type: $optlevel; exit; break;;
esac

case "$realprec" in
   single) precmap="{'real':{'':'float'}}"; shift;;
   double) compflags+=(-fdefault-real-8); precmap="{'real':{'':'double'}}"; shift;;
   *) echo Invalid precision type: $realprec; exit; break;;
esac

if [[ -d $BINDIR ]]; then
   case $buildtype in
   program) rm -f "$BINDIR"/"$NAME";;
   library) rm -f "$BINDIR"/"$NAME".so "$BINDIR"/"$NAME".*.so;;
   esac
else
   mkdir "$BINDIR"
fi

if [[ -d $BUILDIR ]]; then
  if $recompile; then
     rm -f "$BUILDIR"/*
  fi
else
  mkdir -p "$BUILDIR"
fi

objectlist=()
exportlist=()

while IFS= read -r line; do
  eval set -- "$line"
  compile "$1"
  if [[ -n $2 ]]; then
     exportlist+=("$BUILDIR"/"$1")
  fi
done < <(grep -v '^#' "$SRCDIR"/compilelist)

case $buildtype in
program)
   echo Linking program...
   $FORTRAN "${linkflags[@]}" -llapack -o "$BINDIR"/"$NAME" "${objectlist[@]}"
   ;;
library)
#   echo Linking shared library...
#   $FORTRAN -shared -L"$LAPACK_PATH" -llapack -o "$BINDIR"/"$NAME".so "${objectlist[@]}"
   echo Linking python library...
   pushd "$BINDIR"
   echo "$precmap" > .f2py_f2cmap
   $F2PY -h "$BUILDIR"/"$NAME".pyf --overwrite-signature -m "$NAME" "${exportlist[@]}" --quiet
   $F2PY -c "$BUILDIR"/"$NAME".pyf -I"$BUILDIR" -L"$LAPACK_PATH" -llapack "${objectlist[@]}" --quiet
   popd
esac
