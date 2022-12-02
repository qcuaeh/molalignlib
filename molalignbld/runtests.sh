#!/bin/bash
shopt -s nullglob

buildir=$1
testdir=$2
shift 2

for file in "$testdir"/*.out; do
   name=$(basename "$file")
   echo -n Running test ${name%.out}.xyz...
   if diff -bB "$file" <("./$buildir/molalign" "$testdir/${name%.out}.xyz" "$@") > /dev/null; then
       echo \ pass
   else
       echo \ failed
   fi
done

for file in aligned_*.xyz; do
   rm "$file"
done
