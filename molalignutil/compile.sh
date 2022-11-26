#!/bin/bash -e
shopt -s nullglob

split() {
   IFS=$2 read -r -a "$1" <<< "${!1}"
}

all=false
pic=false
debug=false
real=8

split base_flags \ 
split optim_flags \ 
split debug_flags \ 
split double_flags \ 

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
   exit 1
fi

if test ! -d "$1"; then
   echo Error: $1 does not exist
   exit 1
fi

srcdir=$(cd "$1"; pwd)

if test ! -e "$2"; then
   mkdir "$2"
elif test ! -d "$2"; then
   echo Error: $2 does exist but is not a directory
   exit 1
fi

pushd "$2" >/dev/null

flags=("${base_flags[@]}")

if $pic; then
   flags+=(-fPIC)
fi

if $debug; then
   flags+=("${debug_flags[@]}")
else
   flags+=("${optim_flags[@]}")
fi

case "$real" in
4)
   echo '{"real":{"":"float"}}' > .f2py_f2cmap
   shift
   ;;
8)
   flags+=("${double_flags[@]}")
   echo '{"real":{"":"double"}}' > .f2py_f2cmap
   shift
   ;;
*)
   echo Invalid precision type: $real
   exit 1
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
