#!/bin/bash -e
shopt -s nullglob
unalias -a

to_array() {
   IFS=\  read -r -a "$1" <<< "${!1}"
}

build_library() {
   if $full_build; then
      full_build=false
      if test -d build; then
         for file in build/*{.f90,.mod,.o}; do
            rm "$file"
         done
      fi
   fi
   if test ! -d "$libdir"; then
      echo Error: $libdir does not exist
      exit 1
   fi
   comp_flags=("${std_flag[@]}")
   comp_flags+=("${extra_flags[@]}")
   if $pic_build; then
      comp_flags+=("${pic_flag[@]}")
   fi
   if $debug_build; then
      comp_flags+=("${debug_flags[@]}")
   else
      comp_flags+=("${optim_flags[@]}")
   fi
   while IFS= read -r srcfile; do
      objfile=${srcfile%.*}.o
      if ! test -e "$buildir/$objfile" \
      || ! test -e "$buildir/$srcfile" \
      || ! diff -q "$buildir/$srcfile" "$libdir/$srcfile" >/dev/null
      then
#         echo $srcfile added to compile list
         cp -f "$libdir/$srcfile" "$buildir/$srcfile"
         compile_list+=("$srcfile")
      fi
      object_files+=("$objfile")
   done < <(grep -v ^# "$libdir/source_files")
   pushd "$buildir" >/dev/null
   for srcfile in "${compile_list[@]}"; do
      echo Recompiling $srcfile...
      "$F90" "${comp_flags[@]}" -c "$srcfile"
   done
   ar r molalignlib.a "${object_files[@]}"
   popd >/dev/null
}

build_program() {
   if test -z "$1"; then
      echo Error: name is empty
      exit 1
   fi
   echo Building program ${1%.*}...
   cp "$rootdir/$1" "$buildir"
   pushd "$buildir" > /dev/null
   "$F90" "${comp_flags[@]}" "${link_flags[@]}" "$1" molalignlib.a -o "${1%.*}"
   popd > /dev/null
}

run_tests() {
   suffix=$1
   subdir=$2
   shift 2
   executable=$buildir/atomalign
   for file in "$testdir/$subdir"/*.xyz; do
      name=$(basename "$file" .xyz)_$suffix
      echo -n "Running test $subdir/$name... "
      if $write_test; then
         "$executable" -pipe -test -stats -N 5 "$@" < "$file" > "$testdir/$subdir/$name.out" 2>/dev/null
         echo done
      else
         if diff -bB "$testdir/$subdir/$name.out" <("$executable" -pipe -test -stats -N 5 "$@" < "$file" 2>/dev/null); then
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
libdir=$rootdir/molalignlib

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

to_array std_flag
to_array extra_flags
to_array pic_flag
to_array optim_flags
to_array debug_flags
to_array link_flags
to_array f2py_flags

# Read arguments and set options

build_flag=true
test_flag=true
write_test=false
full_build=true
debug_build=false
pic_build=false

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

if $build_flag; then
   # Build program
   build_library
   build_program atomalign.f90
   build_program molalign.f90
fi

if $test_flag; then
   # Run tests
   run_tests prune17 jcim.2c01187/0.05 -remap -prune rd -tol 0.17
#   run_tests bondbiasmna MOBH35-shuffled -remap -bond -bias mna
#   run_tests bondbiasmnaback MOBH35-shuffled -remap -bond -bias mna -back
#   run_tests bondbiasmnabackreac MOBH35-shuffled -remap -bond -bias mna -back -reac
fi
