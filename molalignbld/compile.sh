#!/bin/bash -e
shopt -s nullglob

split() {
   IFS=$2 read -r -a "$1" <<< "${!1}"
}

pic=false
debug=false

split base_flags \ 
split optim_flags \ 
split debug_flags \ 
split double_flags \ 

options=$(getopt -o '' -al pic,debug -- "$@") || exit
eval set -- "$options"

while true; do
   case "$1" in
   --pic) pic=true; shift ;;
   --debug) debug=true; shift ;;
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

while IFS= read -r srcfile; do
   objfile=${srcfile%.f*}.o
   if ! test -e "$objfile" \
   || ! test -e "$srcfile" \
   || ! diff -q "$srcfile" "$srcdir/$srcfile" >/dev/null
   then
      echo Compiling ${srcfile%.f*}...
      if test "$srcdir" != "$PWD"; then
         cp -f "$srcdir/$srcfile" "$srcfile"
      fi
      "$F90" "${flags[@]}" -c "$srcfile" -o "$objfile"
   fi
done < <(grep -v '^#' "$srcdir/fortran_files")

if test -f "$srcdir/f2py_files"; then
   while IFS= read -r srcfile; do
      cp -f "$srcdir/$srcfile" "${srcfile%.f*}.f2py"
   done < <(grep -v '^#' "$srcdir/f2py_files")
fi

popd >/dev/null
