#!/bin/sh
set -e

all_args="$@"
mpi="false"
travis="false"
cmake_args=""
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

    case "$1" in -D*)
                     cmake_args="$cmake_args $1"
                     shift
                     continue
                     ;;
    esac

    # parse pairs
    eval $( echo $1 | sed 's/-//g' | tr -d '\012')=$2
    shift
    shift
done

if [ "$mpi" = "true" ] && [ "$travis" = "false" ]; then
    BUILD_DIR="build-mpi"
else
    BUILD_DIR="build"
fi

# This is bad because it overrides any choices the user may have chosen.
CC=gcc CXX=g++

if [ "$travis" = "true" ]; then
    CC=${C_COMPILER}
    CXX=${CXX_COMPILER}
    all_args="-travis true -mpi ${USE_MPI} -help true -exec_name rb"
fi

export CC
export CXX
echo "CC=${CC}  CXX=${CXX}"

if [ "$1" = "clean" ]
then
	rm -rf ${BUILD_DIR}
else
if [ ! -d ${BUILD_DIR} ]; then
	mkdir ${BUILD_DIR}
fi

    #################
    # generate git version number
    ./generate_version_number.sh
    if [ -e ../../src/revlanguage/utils/GitVersion.cpp ] ; then
        cp ../../src/revlanguage/utils/GitVersion.cpp GitVersion_backup.cpp
    fi
    mv GitVersion.cpp ../../src/revlanguage/utils/

	./regenerate.sh ${all_args}
	cd ${BUILD_DIR}
        echo "Running 'cmake . $cmake_args' in $(pwd)"
	cmake . $cmake_args
        echo
        echo "Running 'make -j4' in $(pwd)"
	make -j 4
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
