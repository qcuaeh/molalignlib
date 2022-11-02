#!/bin/bash
# running and comparing tests with outputs stored in tests/outputs/*/*.out

indir=tests/r01
outdir=tests/outputs/r01
opts="-sort -bias -tol 0.35 -test -rec 10"

echo
echo "Inputs are in directory $indir"
echo "Reference outputs are in directory $outdir"
echo "Tests will be run with options $opts"
echo

for file in $( ls $outdir/*.out ); do
    name=$( basename $file )
    echo -n "Running test ${name%.out}.xyz... "
    if diff -bB $file <(bin/molalign tests/r01/${name%.out}.xyz $opts) > /dev/null; then
        echo "ok"
    else
        echo "Failed"
    fi
done

rm aligned_*.xyz
echo "Done"
echo

