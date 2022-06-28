# Compile source file
compile () {
   sourcefile=$SRCDIR/$1
   buildfile=$OBJDIR/$1
   objectfile=$OBJDIR/${1%.*}.o
   objectlist+=("$objectfile")
   if ! test -e "$objectfile" || ! test -e "$buildfile" || ! diff -q "$sourcefile" "$buildfile" > /dev/null; then
      rm -f "$objectfile"
      cp -p "$sourcefile" "$buildfile"
      echo Compiling "$1"
      $FORTRAN $libflags $optflags -c "$sourcefile" -o "$objectfile" -I "$OBJDIR" -J "$OBJDIR" || exit
   fi
}

shopt -s nullglob

options=$(getopt -a -o '' -l lib,slow,fast,debug -- "$@")
if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi
eval set -- "$options"

optlevel=fast
buildtype=executable

while true; do
   case "$1" in
      --slow) optlevel=slow; shift;;
      --fast) optlevel=fast; shift;;
      --debug) optlevel=debug; shift;;
      --lib) buildtype=library; shift;;
      --) shift; break ;;
   esac
done

case $buildtype in
   executable) libflags='';;
   library) libflags='-fPIC';;
   *) echo Invalid build type: $buildtype; exit;;
esac

case $optlevel in
   slow) optflags='-O3';;
   fast) optflags='-O3 -ffast-math';;
   debug) optflags='-g -fbounds-check -fbacktrace -ffpe-trap=zero,invalid,overflow -O0 -Wall';;
   *) echo Invalid optimization level: $optlevel; exit;;
esac

NAME=ralign
ROOTDIR=$(cd -- "$(dirname "$0")" && pwd)
SRCDIR=$ROOTDIR/source
BINDIR=$ROOTDIR/bin
BUILDIR=$ROOTDIR/_build_dir
OBJDIR=$BUILDIR/$buildtype/$optlevel

if [[ -d $BINDIR ]]; then
   case $buildtype in
      executable) rm -f "$BINDIR"/"$NAME";;
      library) rm -f "$BINDIR"/"$NAME".so "$BINDIR"/"$NAME".*.so;;
   esac
else
   mkdir "$BINDIR"
fi

if [[ -d $OBJDIR ]]; then
  if $RECOMPILE; then
     rm -f "$OBJDIR"/*
  fi
else
  mkdir -p "$OBJDIR"
fi

objectlist=()
exportlist=()

while IFS= read -r line; do
  eval set -- "$line"
  compile "$1"
  if [[ -n $2 ]]; then
     exportlist+=("$OBJDIR"/"$1")
  fi
done < <(grep -v '^#' "$SRCDIR"/compilelist)

case $buildtype in
executable)
   echo Linking executable...
   $FORTRAN -L "$LAPACK" -llapack -o "$BINDIR"/"$NAME" "${objectlist[@]}"
   ;;
library)
#   echo Linking shared library...
#   $FORTRAN -shared -L "$LAPACK" -llapack -o "$BINDIR"/"$NAME".so "${objectlist[@]}"
   echo Linking python library...
   pushd "$BINDIR"
   $F2PY -h "$OBJDIR"/"$NAME".pyf --overwrite-signature -m "$NAME" "${exportlist[@]}" --quiet
   $F2PY -c "$OBJDIR"/"$NAME".pyf -I"$OBJDIR" -L"$LAPACK" -llapack "${objectlist[@]}" --quiet
   popd
esac
