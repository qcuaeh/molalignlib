#!/bin/bash -e
shopt -s nullglob
unalias -a

to_array() {
   IFS=\  read -r -a "$1" <<< "${!1}"
}

build_library() {
   srcdir=$PWD/molalignlib
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
   while read -r filename; do
      srcfile=$srcdir/$filename
      bldfile=$blddir/$filename
      objfile=$blddir/${filename%.*}.o
      if ! test -e "$objfile" \
      || ! test -e "$bldfile" \
      || ! diff -q "$bldfile" "$srcfile" >/dev/null
      then
         cp "$srcfile" "$bldfile"
         echo Compiling $filename...
         "$F90" "${comp_flags[@]}" -J "$blddir" -c "$bldfile" -o "$objfile"
      fi
      object_files+=("$objfile")
   done < <(grep -v ^\# "$srcdir/source_files")
   ar r "$blddir/molalignlib.a" "${object_files[@]}"
}

build_programs() {
   srcdir=$PWD/programs
   molalignlib=$blddir/molalignlib.a
   while read -r filename progname; do
      srcfile=$srcdir/$filename
      bldfile=$blddir/$filename
      execfile=$blddir/$progname
      cp "$srcfile" "$bldfile"
      echo Building program $progname...
      "$F90" "${comp_flags[@]}" "${link_flags[@]}" "$bldfile" "$molalignlib" -o "$execfile"
   done < <(grep -v ^\# "$srcdir/program_files")
}

run_tests() {
   suffix=$1
   subdir=$2
   shift 2
   testdir=$PWD/tests
   executable=$blddir/rmsd-cluster
   for file in "$testdir/$subdir"/*.xyz; do
      name=$(basename "$file" .xyz)_$suffix
      echo -n "Running test $subdir/$name... "
      if $write_test; then
         "$executable" "$@" < "$file" > "$testdir/$subdir/$name.out" 2>/dev/null
         echo done
      else
         if diff -bB "$testdir/$subdir/$name.out" <("$executable" "$@" < "$file" 2>/dev/null); then
            echo ok
         else
            echo failed
         fi
      fi
   done
}

if test ! -e ./build.env; then
   echo Error: build.env does not exist
   exit 1
elif test ! -f ./build.env; then
   echo Error: build.env does exist but is not a file
   exit 1
fi

blddir=$PWD/build

if test ! -e "$blddir"; then
   mkdir "$blddir"
elif test ! -d "$blddir"; then
   echo Error: $blddir does exist but is not a directory
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

# Build static library and programs
if $build_flag; then
   build_library
   build_programs
fi

# Run tests
if $test_flag; then
   run_tests prune17 jcim.2c01187/0.05 -pipe xyz -test -stats -N 5 -align -remap -prune rd -tol 0.17
#   run_tests bondbiasmna MOBH35-shuffled -remap -bond -bias mna
#   run_tests bondbiasmnaback MOBH35-shuffled -remap -bond -bias mna -back
#   run_tests bondbiasmnabackreac MOBH35-shuffled -remap -bond -bias mna -back -reac
fi
