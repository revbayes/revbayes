#!/bin/bash


if [ -z "$1" ] ; then
    printf "Please supply the full path to rb as first argument.\n\n"
    printf "Examples:\n"
    printf '  ./run_integration_tests.sh "$(readlink -f ../projects/cmake/build/rb)"\n'
    printf '  ./run_integration_tests.sh "$PWD/../projects/cmake/build/rb"\n'
    printf '  ./run_integration_tests.sh  -mpi true "$PWD/../projects/cmake/build-mpi/rb-mpi"\n'
#    printf '  ./run_integration_tests.sh mpirun -np 4 "$(readlink -f ../projects/cmake/rb)"\n'
    exit 101
fi

shopt -s nullglob  # make *.Rev expand to the empty string if there are no matches.

RED="\033[31;1m" #bold red
GREEN="\033[32;1m" #bold green
BLUE="\033[34;1m" # bold blue
BOLD="\033[1m"
CLEAR="\033[0m"
UNDERLINE="\033[4m"

BLUE2=$'\033[34;1m'
CLEAR2=$'\033[0m'


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
    echo "RevBayes command '${rb_exec}' seems not to work!"
    exit 102
fi

tests=()
status=()

if [ ! -d "revbayes.github.io" ] ; then
    echo "No revbayes.github.io directory, cloning it"
    git clone https://github.com/revbayes/revbayes.github.io.git
fi

# Run the tutorial tests using the script from the website
(
    cd revbayes.github.io/tutorials
    ./run_tutorial_tests.sh ${rb_exec}
)

for t in test_*; do
    testname=`echo $t | cut -d _ -f 2-`

    if [ -d $t ]; then
        tests+=($testname)
    else
        continue
    fi

    printf "\n${BOLD}#### Running test: ${CLEAR}${UNDERLINE}$testname${CLEAR}\n"
    cd $t

    rm -rf output data
    ln -s ../data data

    test_result=0
    # run the test scripts
    for f in scripts/*.[Rr]ev ; do
        printf "   ${f}: "
        mkdir -p output
        tmp0=${f#scripts/}
        script_name=${tmp0%.[Rr]ev}
        ${rb_exec} -b $f &> output/${script_name}.errout # print output so we can see any error messages
        script_result="$?"

        expected_exit=0
        if [ -f "exit_expected/${script_name}" ] ; then
            expected_exit="$(cat exit_expected/${script_name})"
        fi
        if [ "${script_result}" = 139 ]; then
            script_result="SEGFAULT"
        elif [ "${script_result}" = 134 ]; then
            script_result="Aborted"
        elif [ "${script_result}" != "${expected_exit}" ]; then
            script_result="error ${script_result}"
        fi

        if [ "${script_result}" != "${expected_exit}" ] ; then
            script_result="${f}: ${script_result}"
            echo
            tail -n 5 "output/${script_name}.errout" | sed "s/^/       ${BLUE2}|${CLEAR2}  /g"
            printf "\n   ${RED}FAIL${CLEAR}: ${script_result}\n"
            echo
            if [ "${test_result}" = 0 ] ; then
                test_result="\t${script_result}\n"
            else
                test_result="${test_result}\t${script_result}\n"
            fi
        else
            printf "${GREEN}done${CLEAR}.\n"
        fi

    done

    # store the exit status
    status+=("${test_result}")

    rm data

    cd $OLDPWD

done

printf "\n\n${BOLD}#### Checking output from tests... ${CLEAR}\n"
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

        # some special handling for the *.errout files
        for f in scripts/*.[Rr]ev ; do
            tmp0=${f#scripts/}
            script_name=${tmp0%.[Rr]ev}

            # Delete all before the 1st occurrence of the string '   Processing file' (inclusive)
            # Use a temporary intermediate file to make this work w/ both GNU and BSD sed
            sed '1,/   Processing file/d' output/${script_name}.errout > output/${script_name}.errout.tmp
            mv output/${script_name}.errout.tmp output/${script_name}.errout

            # Also delete the final line of failing tests, which reprints the path to the script
            # that differs between Windows and Unix (has no effect if the line is absent)
            sed '/   Error:\tProblem processing/d' output/${script_name}.errout > output/${script_name}.errout.tmp
            mv output/${script_name}.errout.tmp output/${script_name}.errout

            # Account for OS-specific differences in path separators
            if [ "$windows" = "true" ]; then
                sed 's/\\/\//g' output/${script_name}.errout > output/${script_name}.errout.tmp
                mv output/${script_name}.errout.tmp output/${script_name}.errout
            fi
        done

        for expected in $(find output_expected -type f 2>/dev/null); do
            filename=${expected#output_expected/}
            output=output/${filename}
            if [ ! -e "${output}" ]; then
                errs+=("missing:  ${filename}")
            elif ! diff "${output}" "${expected}" > /dev/null; then
                errs+=("mismatch: ${filename}")
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
            printf "${BOLD}>>>> Test ${MAGENTA}failed${CLEAR}${BOLD}: ${CLEAR}${UNDERLINE}$t${CLEAR} (expected)\n"
        else
            ((failed++))
            printf "${BOLD}>>>> Test ${RED}failed${CLEAR}${BOLD}: ${CLEAR}${UNDERLINE}$t${CLEAR}\n"
        fi
        for errmsg in "${errs[@]}"; do
            printf "\t$errmsg\n"
        done
    else
        ((pass++))
        printf "${BOLD}#### Test ${GREEN}passed${CLEAR}${BOLD}: ${CLEAR}${UNDERLINE}$t${CLEAR}\n"
    fi

    ((i++))
done


if [ $failed -gt 0 ]; then
    printf "\n\n${BOLD}#### ${MAGENTA}Warning!${CLEAR}${BOLD} unexpected failures: $failed   expected failures: $xfailed   total tests: $i${CLEAR}\n\n"
    exit 113
else
    printf "\n\n${BOLD}#### ${GREEN}Success!${CLEAR}${BOLD} unexpected failures: $failed   expected failures: $xfailed   total tests: $i${CLEAR}\n\n"
    printf "\n\n${BOLD}#### All tests passed.${CLEAR}\n\n"
fi
