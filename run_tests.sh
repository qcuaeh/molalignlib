#!/bin/bash
# running and comparing tests with outputs stored in tests/outputs

run_test() {
    echo Inputs are in directory $inputdir
    echo Reference outputs are in directory $outputdir
    echo Tests will be run with options ${options[@]}
    echo

    for file in "$outputdir"/*.out; do
        name=$(basename "$file")
        echo -n "Running test ${name%.out}.xyz... "
        if diff -bB "$file" <(./molalign "$inputdir/${name%.out}.xyz" "${options[@]}") > /dev/null; then
            echo ok
        else
            echo Failed
        fi
    done

    rm aligned_*.xyz
    echo "Done"
    echo
}


cd "$(dirname "$0")"

if [[ -f build.env ]]; then
   . build.env
else
   echo Error: build.env does not exist or is not a file
fi

inputdir=tests/r005
outputdir=tests/outputs/r005
options=(-test -rec 10 -sort -bias -tol 0.17)
run_test

inputdir=tests/r01
outputdir=tests/outputs/r01
options=(-test -rec 10 -sort -bias -tol 0.35)
run_test

inputdir=tests/r02
outputdir=tests/outputs/r02
options=(-test -rec 10 -sort -bias -tol 0.69)
#run_test

