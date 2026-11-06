#!/bin/bash -e
shopt -s nullglob
unalias -a

run_tests() {
   testdir=$1
   options=(-align -remap -printstats -maxrecs 10 -prune "$1" -coords xyz)
   while read name; do
      echo -n "Running test $testdir/$name... "
      if $write_test; then
         atormsd "${options[@]}" "$testdir/${name}_1.xyz" "$testdir/${name}_2.xyz" > "$testdir/$name.out" 2>&1
         echo done
      else
         if diff -bB "$testdir/$name.out" <(atormsd "${options[@]}" "$testdir/${name}_1.xyz" "$testdir/${name}_2.xyz" 2>&1); then
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

run_tests 0.05
#run_tests 0.1
#run_tests 0.2
