#!/bin/bash -e
shopt -s nullglob
unalias -a

test_binary() {
  if ! test -e "$1"; then
    echo $1 does not exist
    exit 1
  fi
  if ! test -x "$1"; then
    echo $1 is not executable
    exit 1
  fi
}

run_tests_jcim2c01187() {
   rmsdbin=./atormsd
   testdir=jcim.2c01187/$1
   options=(-align -remap -stats -records 10 -prune "$1" -aligned xyz)
   test_binary "$rmsdbin"
   while read name; do
      echo -n "Running test $testdir/$name... "
      if $write_test; then
         $rmsdbin "${options[@]}" "$testdir/${name}_1.xyz" "$testdir/${name}_2.xyz" > "$testdir/$name.out" 2>&1
         echo done
      else
         if diff -bB "$testdir/$name.out" <($rmsdbin "${options[@]}" "$testdir/${name}_1.xyz" "$testdir/${name}_2.xyz" 2>&1); then
            echo ok
         else
            echo failed
         fi
      fi
   done < $testdir/test_files
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

run_tests_jcim2c01187 0.05
#run_tests_jcim2c01187 0.1
#run_tests_jcim2c01187 0.2
