shopt -s nullglob

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
LIBDIR=$HERE/lib
BUILDROOT=$HERE/_build_dir

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

BUILDIR=$BUILDROOT/$buildtype/$realprec/$optlevel

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
   single) shift;;
   double) compflags+=(-fdefault-real-8); shift;;
   *) echo Invalid precision type: $realprec; exit; break;;
esac

case $buildtype in
program)
   if [[ -d $BINDIR ]]; then
      rm -f "$BINDIR"/"$NAME"
   else
      mkdir "$BINDIR"
   fi
   ;;
library)
   if [[ -d $LIBDIR ]]; then
      rm -f "$LIBDIR"/"$NAME".so "$LIBDIR"/"$NAME".*.so
   else
      mkdir "$LIBDIR"
   fi
   ;;
esac

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
   $FORTRAN "${libpathlist[@]}" -llapack -o "$BINDIR"/"$NAME" "${objectlist[@]}"
   ;;
library)
   echo Linking libraries...
#   $FORTRAN -shared "${libpathlist[@]}" -llapack -o "$LIBDIR"/"$NAME".so "${objectlist[@]}"
   pushd "$LIBDIR"
   $F2PY -h "$BUILDIR"/"$NAME".pyf --overwrite-signature -m "$NAME" "${exportlist[@]}" --quiet
   $F2PY -c "$BUILDIR"/"$NAME".pyf -I"$BUILDIR" "${libpathlist[@]}" -llapack "${objectlist[@]}" --quiet
   popd
esac
