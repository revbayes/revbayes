#!/bin/sh
set -e

all_args="$@"

#################
# command line options
# set default values
debug="false"
travis="false"
mpi="false"
gentoo="false"
help="false"
jupyter="false"
boost_root=""
boost_lib=""
beagle="false"
static_boost="false"
j=4

cmake_args=""
# parse command line arguments
while echo $1 | grep ^- > /dev/null; do
    # intercept help while parsing "-key value" pairs
    if [ "$1" = "--help" ] || [ "$1" = "-h" ]
    then
        echo 'Command line options are:
-h                              : print this help and exit.
-debug          <true|false>    : set to true to build in debug mode. Defaults to false.
-mpi            <true|false>    : set to true if you want to build the MPI version. Defaults to false.
-cmd            <true|false>    : set to true if you want to build RevStudio with GTK2+. Defaults to false.
-jupyter        <true|false>    : set to true if you want to build the jupyter version. Defaults to false.
        (
            cd ../../src
            echo "Generating help database"
            perl ../help/md2help.pl ../help/md/*.md > core/help/RbHelpDatabase.cpp
        )
    fi

    echo "Running './regenerate.sh $(pwd)/$BUILD_DIR"
    ./regenerate.sh $(pwd)/$BUILD_DIR
    cd ${BUILD_DIR}
    echo
    echo
    echo "Running 'cmake ../../../src $cmake_args' in $(pwd)"
    cmake ../../../src $cmake_args
    echo
    echo "Running 'make -j $j' in $(pwd)"
    make -j $j
    cd ..

    if [ -e  GitVersion_backup.cpp ] ; then
        cp GitVersion_backup.cpp ../../src/revlanguage/utils/GitVersion.cpp
        rm GitVersion_backup.cpp
    fi
fi

# Run tests
if [ "$travis" = "true" ] && [ "${TRAVIS_BUILD_STAGE_NAME}" = "Test" ]
then
  cd ../..
  echo "\"Hello World\"" | projects/cmake/rb
  cd tests
  ./run_integration_tests.sh -mpi ${USE_MPI} ${TRAVIS_BUILD_DIR}/projects/cmake/rb
  # Run testiphy
  export PATH=${TRAVIS_BUILD_DIR}/projects/cmake:$PATH
  cd
  git clone https://gitlab.com/testiphy/testiphy.git
  cd testiphy
  ./testiphy rb
fi
