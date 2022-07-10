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
      $FORTRAN $libflags $optflags -c "$sourcefile" -o "$objectfile" -I "$BUILDIR" -J "$BUILDIR" || exit
   fi
}

shopt -s nullglob

options=$(getopt -a -o '' -l slow,fast,debug,recompile -- "$@")
if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi
eval set -- "$options"

optlevel=fast
recompile=false

while true; do
   case "$1" in
      --slow) optlevel=slow; shift;;
      --fast) optlevel=fast; shift;;
      --debug) optlevel=debug; shift;;
      --recompile) recompile=true; shift;;
      --) shift; break ;;
   esac
done

case $BUILDTYPE in
   program) libflags='';;
   library) libflags='-fPIC';;
   *) echo Invalid build type: $BUILDTYPE; exit;;
esac

case $optlevel in
   slow) optflags='-O0';;
   fast) optflags='-O3 -ffast-math';;
   debug) optflags='-O0 -g -fbounds-check -fbacktrace -Wall -ffpe-trap=zero,invalid,overflow';;
   *) echo Invalid optimization level: $optlevel; exit;;
esac

NAME=ralign
HERE=$(cd -- "$(dirname "$0")" && pwd)
SRCDIR=$HERE/source
BINDIR=$HERE/bin
BUILDROOT=$HERE/_build_dir
BUILDIR=$BUILDROOT/$BUILDTYPE/$optlevel

if [[ -d $BINDIR ]]; then
   case $BUILDTYPE in
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

case $BUILDTYPE in
program)
   echo Linking program...
   $FORTRAN -L "$LAPACK" -llapack -o "$BINDIR"/"$NAME" "${objectlist[@]}"
   ;;
library)
#   echo Linking shared library...
#   $FORTRAN -shared -L "$LAPACK" -llapack -o "$BINDIR"/"$NAME".so "${objectlist[@]}"
   echo Linking python library...
   pushd "$BINDIR"
   $F2PY -h "$BUILDIR"/"$NAME".pyf --overwrite-signature -m "$NAME" "${exportlist[@]}" --quiet
   $F2PY -c "$BUILDIR"/"$NAME".pyf -I"$BUILDIR" -L"$LAPACK" -llapack "${objectlist[@]}" --quiet
   popd
esac
