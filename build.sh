NAME=ralign
F2PY=f2py3.4
FORTRAN=gfortran
LIBPATH=/usr/lib64/atlas

# Compile source file
compile () {
   sourcefile=$SRCDIR/$1
   buildfile=$OBJDIR/$1
   objectfile=$OBJDIR/${1%.*}.o
   objectlist+=($objectfile)
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
buildtype=program

while true; do
   case "$1" in
      --slow) optlevel=slow; shift;;
      --fast) optlevel=fast; shift;;
      --debug) optlevel=debug; shift;;
      --lib) buildtype=library; shift;;
      --) shift; break ;;
   esac
done

case "$buildtype" in
   program) libflags='';;
   library) libflags='-fPIC';;
   *) echo Invalid build type: $buildtype; exit;;
esac

case "$optlevel" in
   slow) optflags='-O3';;
   fast) optflags='-O3 -ffast-math';;
   debug) optflags='-g -fbounds-check -fbacktrace -ffpe-trap=zero,invalid,overflow -O0 -Wall';;
   *) echo Invalid optimization level: $optlevel; exit;;
esac

BINDIR=$(dirname "$0")/bin
SRCDIR=$(dirname "$0")/fortran
SRCLIST=$(dirname "$0")/srclist
BUILDIR=$(dirname "$0")/_build_dir
OBJDIR=$BUILDIR/$buildtype/$optlevel

if [[ -d $BINDIR ]]; then
   for file in "$BINDIR"/*; do
      rm "$file"
   done
else
   mkdir "$BINDIR"
fi

if [[ -d $OBJDIR ]]; then
  if $RECOMPILE; then
     rm -r "$OBJDIR"/*
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
     exportlist+=($OBJDIR/$1)
  fi
done < <(grep -v '^#' "$SRCLIST")

case "$buildtype" in
   program)
      echo Linking program...
      $FORTRAN -L$LIBPATH -llapack -o "$BINDIR/$NAME" "${objectlist[@]}"
      ;;
   library)
      echo Linking library...
      $F2PY --quiet --overwrite-signature -m $NAME -h $NAME.pyf "${exportlist[@]}"
      $F2PY --quiet --build-dir $OBJDIR -I$OBJDIR -L$LIBPATH -llapack -c $NAME.pyf "${objectlist[@]}"
      ;;
esac
