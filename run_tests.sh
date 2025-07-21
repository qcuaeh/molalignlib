#!/bin/bash -e
shopt -s nullglob
unalias -a

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

run_tests prune17 jcim.2c01187/0.05 -pipe xyz -test -stats -N 5 -align -remap -prune rd -tol 0.17
