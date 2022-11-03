#!/bin/bash
# running and comparing tests with outputs stored in tests/outputs/*/*.out

run_test() {
    echo "Inputs are in directory $indir"
    echo "Reference outputs are in directory $outdir"
    echo "Tests will be run with options ${opts[@]}"
    echo

    for file in $( ls $outdir/*.out ); do
        name=$( basename $file )
        echo -n "Running test ${name%.out}.xyz... "
        if diff -bB $file <(bin/molalign "$indir/${name%.out}.xyz" "${opts[@]}") > /dev/null; then
            echo "ok"
        else
            echo "Failed"
        fi
    done

    rm aligned_*.xyz
    echo "Done"
    echo
}

indir=tests/r005
outdir=tests/outputs/r005
opts=(-test -rec 10 -sort -bias -tol 0.17)
run_test

indir=tests/r01
outdir=tests/outputs/r01
opts=(-test -rec 10 -sort -bias -tol 0.35)
run_test

indir=tests/r02
outdir=tests/outputs/r02
opts=(-test -rec 10 -sort -bias -tol 0.69)
#run_test

