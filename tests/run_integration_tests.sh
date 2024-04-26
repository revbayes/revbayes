#!/bin/bash

git submodule update --init --recursive

if [ -z "$1" ] ; then
    printf "Please supply the full path to rb as first argument.\n\n"
    printf "Examples:\n"
    printf '  ./run_integration_tests.sh "$(readlink -f ../projects/cmake/rb)"\n'
    printf '  ./run_integration_tests.sh "$PWD/../projects/cmake/rb"\n'
    printf '  ./run_integration_tests.sh  -mpi true "$PWD/../projects/cmake/rb"\n'
#    printf '  ./run_integration_tests.sh mpirun -np 4 "$(readlink -f ../projects/cmake/rb)"\n'
    exit 101
fi


mpi="false"

# parse command line arguments
while echo $1 | grep ^- > /dev/null; do
    # intercept help while parsing "-key value" pairs
    if [ "$1" = "--help" ] || [ "$1" = "-h" ]
    then
        echo '
Command line options are:
-h                              : print this help and exit.
'
        exit
    fi

    # parse pairs
    eval $( echo $1 | sed 's/-//g' | tr -d '\012')=$2
    shift
    shift
done

if [ $mpi = "true" ]; then
#    rb_exec="mpirun --oversubscribe -np 4 $@"
#    rb_exec="mpirun -np 4 $@"
    rb_exec="$@"
else
    rb_exec="$@"
fi

export rb_exec


if ! ${rb_exec} --help > /dev/null 2>&1 ; then
    echo "RevBayes command '${rb_exec}' seems not to work!\n"
    exit 102
fi

tests=()
status=()

for t in revbayes.github.io/tutorials/*/test.sh; do
    testname=`echo $t | cut -d '/' -f 2-3`
    dirname=`echo $t | cut -d '/' -f 1-3`
    
    cd $dirname

    tests+=($testname)

    printf "\n\n#### Running test: $testname\n\n"
    sh test.sh
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

    status+=("$res")

    cd -
done

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
        ${rb_exec} -b $f # print output so we can see any error messages
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

    if [ -d test_$t ]; then
        cd test_$t

        # check if output matches expected output
        errs=()
        exp_out_dir="output_expected"
        # Sebastian: For now we only use single cores until we fix Travis mpirun
    #    if [ "${mpi}" = "true" ]; then
    #        if [ -d "output_expected_mpi" ]; then
    #            exp_out_dir="output_expected_mpi"
    #        fi
    #    fi
        if [ "$windows" = "true" ]; then
            find output -type f -exec dos2unix {} \;
            find output -type f -exec sed -i 's/e-00/e-0/g' {} \;
            find output -type f -exec sed -i 's/e+00/e+0/g' {} \;
        fi
        for f in $(ls ${exp_out_dir}); do
            if [ ! -e output/$f ]; then
                errs+=("missing:  $f")
            elif ! diff output/$f ${exp_out_dir}/$f > /dev/null; then
                errs+=("mismatch: $f")
            fi
        done

        cd ..
    fi
    
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

    ((i++))
done


if [ $failed -gt 0 ]; then
    printf "\n\n#### Warning! unexpected failures: $failed   expected failures: $xfailed   total tests: $i\n\n"
    exit 113
else
    printf "\n\n#### Success! unexpected failures: $failed   expected failures: $xfailed   total tests: $i\n\n"
    printf "\n\n#### All tests passed.\n\n"
fi
