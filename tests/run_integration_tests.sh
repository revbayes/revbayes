#!/bin/bash

if [ -z "$1" ] ; then
    printf "Please supply the full path to rb as first argument.\n\n"
    printf "Examples:\n"
    printf '  ./run_integration_tests.sh "$(readlink -f ../projects/cmake/rb)"\n'
    printf '  ./run_integration_tests.sh "$PWD/../projects/cmake/rb"\n'
    printf '  ./run_integration_tests.sh mpirun -np 4 "$(readlink -f ../projects/cmake/rb)"\n'
    exit 101
fi
    

if ! "$@" --help > /dev/null 2>&1 ; then
    echo "RevBayes command '$@' seems not to work!\n"
    exit 102
fi

tests=()
status=()

for t in test_*; do
    testname=`echo $t | cut -d _ -f 2-`

    if [ -d $t ]; then
        tests+=($testname)
    else
        continue
    fi

    printf "\n\n#### Running test: $testname\n\n"
    cd $t

    rm -rf output data
    ln -s ../data data

    res=0
    # run the test scripts
    for f in scripts/*; do
        "$@" -b $f # print output so we can see any error messages
        res="$?"
        if [ $res = 1 ]; then
            res="error: $f"
            break
        elif [ $res = 139 ]; then
            res="segfault: $f"
            break
        elif [ $res != 0 ]; then
            res="error $res: $f"
            break
        fi
    done

    # store the exit status
    status+=("$res")

    rm data

    cd ..
done

printf "\n\n#### Checking output from tests... \n"
xfailed=0
failed=0
pass=0
i=0
while [  $i -lt ${#tests[@]} ]; do
    t=${tests[$i]}
    cd test_$t

    # check if output matches expected output
    errs=()
    for f in $(ls output_expected); do
        if [ ! -e output/$f ]; then
            errs+=("missing:  $f")
        elif ! diff output/$f output_expected/$f > /dev/null; then
            errs+=("mismatch: $f")
        fi
    done

    # check if a script exited with an error
    if [ "${status[$i]}" != 0 ]; then
        errs=("${status[$i]}")
    fi

    # failure if we have an error message
    if [ ${#errs[@]} -gt 0 ]; then
        if [ -f XFAIL ] ; then
            ((xfailed++))
            printf ">>>> Test failed: $t (expected)\n"
        else
            ((failed++))
            printf ">>>> Test failed: $t\n"
        fi
        for errmsg in "${errs[@]}"; do
            printf "\t$errmsg\n"
        done
    else
        ((pass++))
        printf "#### Test passed: $t\n"
    fi

    cd ..

    ((i++))
done


if [ $failed -gt 0 ]; then
    printf "\n\n#### Warning! unexpected failures: $failed   expected failures: $xfailed   total tests: $i\n\n"
    exit 113
else
    printf "\n\n#### Success! unexpected failures: $failed   expected failures: $xfailed   total tests: $i\n\n"
    printf "\n\n#### All tests passed.\n\n"
fi
