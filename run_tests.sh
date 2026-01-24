#!/bin/bash -e
shopt -s nullglob
unalias -a

run_tests_jcim2c01187() {
   subtestdir=$testdir/jcim.2c01187/$1
   options=(-align -remap -test -stats -records 10 -prune "$1")
   rmsdexec=./build/atormsd
   while read name; do
      echo -n "Running test $subtestdir/$name... "
      if $write_test; then
         $rmsdexec "${options[@]}" "$subtestdir/${name}_1.xyz" "$subtestdir/${name}_2.xyz" > "$subtestdir/$name.out" 2>&1
         echo done
      else
         if diff -bB "$subtestdir/$name.out" <($rmsdexec "${options[@]}" "$subtestdir/${name}_1.xyz" "$subtestdir/${name}_2.xyz" 2>&1); then
            echo ok
         else
            echo failed
         fi
      fi
   done < $subtestdir/test_files
}

testdir=$PWD/tests
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

run_tests_jcim2c01187 0.05
#run_tests_jcim2c01187 0.1
#run_tests_jcim2c01187 0.2
