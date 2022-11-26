#!/bin/bash
shopt -s nullglob

run_test() {
    echo Test files are in directory $testdir
    echo Tests will be run with options ${options[@]}
    echo

    for file in "$testdir"/*.out; do
        name=$(basename "$file")
        echo -n Running test ${name%.out}.xyz...
        if diff -bB "$file" <(./molalign "$testdir/${name%.out}.xyz" "${options[@]}") > /dev/null; then
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

testdir=r005
options=(-test -rec 10 -sort -fast -tol 0.17)
run_test

testdir=r01
options=(-test -rec 10 -sort -fast -tol 0.35)
run_test

testdir=tests/r02
options=(-test -rec 10 -sort -fast -tol 0.69)
#run_test
