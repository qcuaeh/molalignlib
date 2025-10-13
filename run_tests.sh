#!/bin/bash -e
shopt -s nullglob
unalias -a

run_tests() {
   suffix=$1
   subdir=$2
   shift 2
   for file in "$testdir/$subdir"/*.xyz; do
      name=$(basename "$file" .xyz)_$suffix
      echo -n "Running test $subdir/$name... "
      if $write_test; then
         "$executable" "$@" < "$file" > "$testdir/$subdir/$name.out" 2>&1
         echo done
      else
         if diff -bB "$testdir/$subdir/$name.out" <("$executable" "$@" < "$file" 2>&1); then
            echo ok
         else
            echo failed
         fi
      fi
   done
}

testdir=$PWD/tests
executable=./build/atormsd
write_test=false

while getopts ":w" opt; do
  case $opt in
    w)
      write_test=true
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

run_tests prune17 jcim.2c01187/0.05 -stdin -coords -stats -N 5 -align -remap -prune 0.05
