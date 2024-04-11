#!/bin/bash -e
shopt -s nullglob
unalias -a

to_array() {
   IFS=\  read -r -a "$1" <<< "${!1}"
}

compile() {
   srcdir=$rootdir/$1
   if $full_build; then
      full_build=false
      if test -d build; then
         for file in build/*{.f90,.mod,.o}; do
            rm "$file"
         done
      fi
   fi
   if test ! -d "$srcdir"; then
      echo Error: $srcdir does not exist
      exit 1
   fi
   pushd "$buildir" >/dev/null
   comp_flags=("${base_flags[@]}")
   if $pic_build; then
      comp_flags+=("${pic_flags[@]}")
   fi
   if $debug_build; then
      comp_flags+=("${debug_flags[@]}")
   else
      comp_flags+=("${opt_flags[@]}")
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
         "$F90" "${comp_flags[@]}" -c "$srcfile" -o "$objfile"
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
   if test -z "$1"; then
      echo Error: name is empty
      exit 1
   fi
   echo Linking program...
   pushd "$buildir" > /dev/null
   "$F90" "${link_flags[@]}" "${obj_files[@]}" -o "$1"
   popd > /dev/null
}

make_lib() {
   if test -z "$1"; then
      echo Error: name is empty
      exit 1
   fi
   echo Linking dynamic library...
   pushd "$buildir" >/dev/null
   "$F90" -shared "${link_flags[@]}" "${obj_files[@]}" -o "$1.so"
   popd >/dev/null
}

make_pyext() {
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
   "$F2PY" -c "$1.pyf" "${f2py_flags[@]}" "${obj_files[@]}" --quiet
   popd >/dev/null
}

run_tests() {
   suffix=$1
   subdir=$2
   shift 2
   executable=$buildir/molalign
   for file in "$testdir/$subdir"/*.xyz; do
      name=$(basename "$file" .xyz)_$suffix
      echo -n "Running test $subdir/$name... "
      if $write_test; then
         "$executable" -pipe -test -stats -N 5 "$@" < "$file" 2>&1 > "$testdir/$subdir/$name.out"
         echo done
      else
         if diff -bB "$testdir/$subdir/$name.out" <("$executable" -pipe -test -stats -N 5 "$@" < "$file" 2>&1); then
            echo ok
         else
            echo failed
         fi
      fi
   done
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

to_array base_flags
to_array pic_flags
to_array opt_flags
to_array debug_flags
to_array link_flags
to_array f2py_flags

# Read arguments and set options

build_flag=true
full_build=true
test_flag=true
write_test=false
debug_build=false

while getopts ":bdqtw" opt; do
  case $opt in
    b)
      build_flag=true
      test_flag=false
      ;;
    d)
      debug_build=true
      ;;
    q)
      full_build=false
      ;;
    t)
      build_flag=false
      test_flag=true
      ;;
    w)
      write_test=true
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

if $build_flag; then
   case $target in
   prog)
      # Build program
      pic_build=false
      compile molalignlib
      compile molalign
      make_prog molalign
      ;;
   lib)
      # Build dynamic library
      pic_build=true
      compile molalignlib
      make_lib molalignlib
      ;;
   pyext)
      # Build python extension module
      pic_build=true
      compile molalignlib
      make_pyext molalignlibext
      ;;
   *)
      echo Unknown target $target
   esac
fi

if $test_flag; then
   # Run tests
   run_tests fast17 jcim.2c01187/0.05 -remap -fast -tol 0.17
   run_tests fastbond MOBH35-shuffled -remap -fast -bond
   run_tests fastbondback MOBH35-shuffled -remap -fast -bond -back
   run_tests fastbondbackreac MOBH35-shuffled -remap -fast -bond -back -reac
fi
