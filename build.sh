#!/bin/bash -e
shopt -s nullglob
unalias -a

toarray() {
   IFS=\  read -r -a "$1" <<< "${!1}"
}

clean_build() {
   $fresh_build || return 0
   fresh_build=false
   if test -d build; then
      pushd build >/dev/null
      for file in *.f90 *.mod *.o; do
         rm "$file"
      done
      popd >/dev/null
   fi
}

compile() {
   $build_flag || return 0
   clean_build
   toarray std_flags
   toarray pic_flags
   toarray optim_flags
   toarray debug_flags
   srcdir=$rootdir/$1
   if test ! -d "$srcdir"; then
      echo Error: $srcdir does not exist
      exit 1
   fi
   pushd "$buildir" >/dev/null
   flags=("${std_flags[@]}")
   if $pic_build; then
      flags+=("${pic_flags[@]}")
   fi
   if $debug_build; then
      flags+=("${debug_flags[@]}")
   else
      flags+=("${optim_flags[@]}")
   fi
   while IFS= read -r srcfile; do
      prefix=${srcfile%.*}
      objfile=$prefix.o
      if ! test -e "$objfile" \
      || ! test -e "$srcfile" \
      || ! diff -q "$srcfile" "$srcdir/$srcfile" >/dev/null
      then
         echo Compiling $prefix...
         if test "$srcdir" != "$PWD"; then
            cp -f "$srcdir/$srcfile" "$srcfile"
         fi
         "$F90" "${flags[@]}" -c "$srcfile" -o "$objfile"
      fi
      obj_files+=("$objfile")
   done < <(grep -v ^# "$srcdir/source_files")
   popd >/dev/null
   if test -f "$srcdir/f2py_files"; then
      while IFS= read -r f2pyfile; do
         f2py_files+=("$f2pyfile")
      done < <(grep -v ^# "$srcdir/f2py_files")
   fi
}

make_prog() {
   executable=$buildir/$1
   $build_flag || return 0
   if test -z "$1"; then
      echo Error: name is empty
      exit 1
   fi
   echo Linking program...
   pushd "$buildir" > /dev/null
   "$F90" -o "$1" "${obj_files[@]}" -llapack
   popd > /dev/null
   echo Done
}

make_lib() {
   $build_flag || return 0
   if test -z "$1"; then
      echo Error: name is empty
      exit 1
   fi
   echo Linking dynamic library...
   pushd "$buildir" >/dev/null
   "$F90" -shared -o "$1.so" "${obj_files[@]}" -llapack
   popd >/dev/null
   echo Done
}

make_pyext() {
   $build_flag || return 0
   if test -z "$1"; then
      echo Error: name is empty
      exit 1
   fi
   if ! type "$F2PY" &>/dev/null; then
      echo Error: F2PY executable not found 
      exit 1
   fi
   pushd "$buildir" >/dev/null
   echo Linking extension module...
   "$F2PY" -h "$1.pyf" -m "$1" --overwrite-signature "${f2py_files[@]}" --quiet
   "$F2PY" -c "$1.pyf" --fcompiler=gnu95 --link-lapack "${obj_files[@]}" --quiet
   popd >/dev/null
   echo Done
}

runtests() {
   $test_flag || return 0
   if test -z "$executable"; then
      echo Error: executable is not set
      exit 1
   fi
   echo "Testing $1"
   subdir=$testdir/$1
   shift
   for file in "$subdir"/*.out; do
      name=$(basename "$file")
      printf "Running test %-16s" ${name%.out}
#      "$executable" -stdin xyz -stdout xyz -test -stats "$@" < "$subdir/${name%.out}.xyz" 2>&1 > "$file"
      if diff -bB <("$executable" -stdin xyz -stdout xyz -test -stats "$@" < "$subdir/${name%.out}.xyz" 2>&1) "$file"; then
          echo Passed
      else
          echo Failed
      fi
   done
   echo Done
}

rootdir=$(dirname "$(readlink -e "$0")")

if test ! -e ./build.env; then
   echo Error: build.env does not exist
   exit 1
elif test ! -f ./build.env; then
   echo Error: build.env does exist but is not a file
   exit 1
fi

buildir=$rootdir/build
testdir=$rootdir/tests

if test ! -e "$buildir"; then
   mkdir "$buildir"
elif test ! -d "$buildir"; then
   echo Error: $buildir does exist but is not a directory
   exit 1
fi

# Set environment
while IFS= read -r line; do
   var=${line%%=*}
   value=${line#*=}
   declare -- "$var"="$value"
done < <(grep -v -e^# -e^$ ./build.env)

build_flag=true
fresh_build=true
debug_build=false
test_flag=true

while getopts ":dqt" opt; do
  case $opt in
    d)
      build_flag=true
      debug_build=true
      test_flag=false
      ;;
    q)
      build_flag=true
      fresh_build=false
      test_flag=false
      ;;
    t)
      build_flag=false
      test_flag=true
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

shift $((OPTIND-1))

if test $# -gt 0; then
   target=$1
else
   target=prog
fi

case $target in
prog)
   # Build program
   pic_build=false
   compile molalignlib
   compile molalign
   make_prog molalign
   # Run tests
   runtests jcim.2c01187/0.05 -rec 5 -sort -fast -tol 0.17
   runtests jcim.2c01187/0.1 -rec 5 -sort -fast -tol 0.35
   runtests jcim.2c01187/0.2 -rec 5 -sort -fast -tol 0.69
   ;;
lib)
   # Build dynamic library
   pic_build=true
   compile molalignlib
   make_lib molalignlib
   echo
   ;;
pyext)
   # Build python extension module
   pic_build=true
   compile molalignlib
   make_pyext molalignlibext
   echo
   ;;
*)
   echo Unknown target $target
esac
