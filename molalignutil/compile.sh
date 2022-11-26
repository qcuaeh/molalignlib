#!/bin/bash -e
shopt -s nullglob

if test ! -e build.env; then
   echo Error: build.env does not exist
   exit
elif test ! -f build.env; then
   echo Error: build.env does exist but is not a file
   exit
fi

. build.env

all=false
pic=false
debug=false
real=8

options=$(getopt -o '' -al all,pic,debug,r4,r8 -- "$@") || exit
eval set -- "$options"

while true; do
   case "$1" in
   --all) all=true; shift ;;
   --pic) pic=true; shift ;;
   --debug) debug=true; shift ;;
   --r4) real=4; shift ;;
   --r8) real=8; shift ;;
   --) shift; break ;;
   *) exit
   esac
done

if test $# -ne 2; then
   echo Error: two arguments are required
   exit
fi

if test ! -d "$1"; then
   echo Error: $1 does not exist
   exit
fi

srcdir=$(cd "$1"; pwd)

if test ! -e "$2"; then
   mkdir "$2"
elif test ! -d "$2"; then
   echo Error: $2 does exist but is not a directory
fi

pushd "$2" >/dev/null

flags=("${BASE_FLAGS[@]}")

if $pic; then
   flags+=(-fPIC)
fi

if $debug; then
   flags+=(-O0 "${DEBUG_FLAGS[@]}")
else
   flags+=(-Ofast)
fi

case "$real" in
4)
   echo '{"real":{"":"float"}}' > .f2py_f2cmap
   shift
   ;;
8)
   flags+=("${REAL8_FLAGS[@]}")
   echo '{"real":{"":"double"}}' > .f2py_f2cmap
   shift
   ;;
*)
   echo Invalid precision type: $real
   exit
   ;;
esac

while IFS= read -r srcfile; do
   objfile=${srcfile%.f*}.o
   if $all \
      || ! test -e "$objfile" \
      || ! test -e "$srcfile" \
      || ! diff -q "$srcfile" "$srcdir/$srcfile" >/dev/null
   then
      echo Compiling ${srcfile%.f*}...
      if test "$srcdir" != "$PWD"; then
         cp -f "$srcdir/$srcfile" "$srcfile"
      fi
      "$F90" "${flags[@]}" -c "$srcfile" -o "$objfile"
   fi
done < <(grep -v '^#' "$srcdir/f90_files")

if test -f "$srcdir/f2py_files"; then
   while IFS= read -r srcfile; do
      cp -f "$srcdir/$srcfile" "${srcfile%.f*}.f2py"
   done < <(grep -v '^#' "$srcdir/f2py_files")
fi

popd >/dev/null
