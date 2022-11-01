#!/bin/bash
# running and comparing tests with outputs stored in tests/outputs/*/*.out

opts="-sort -bias -tol 0.35 -test -rec 10"
perts="r005 r01 r02"
outputDir="outputs"

echo
echo "Checking identical results for new runs of tests..."
echo
echo "Running molalign tests using options: $opts"
echo "Reference outputs in dir: $outputDir; Perturbations: $perts"
echo

for pert in $perts; do
    echo "Tests with perturbation level: $pert"
    for file in $( ls $outputDir/$pert/*.out ); do
        name=$( basename $file )
        echo -n "    Processing test ${name%.out}.xyz ... "
        ../bin/molalign ${pert}/${name%.out}.xyz $opts > tmp.out
        diff=$( diff $file tmp.out )
        if [ "$diff" == "" ]; then res="Same result!"; else res="Different result!"; fi
        echo $res
    done
done

rm aligned_*.xyz tmp.out
echo "Done..."

