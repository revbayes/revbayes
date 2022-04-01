#!/bin/sh
set -e

all_args="$@"

#################
# command line options
# set default values
debug="false"
travis="false"
mpi="false"
help="false"
jupyter="false"
boost_root=""
boost_lib=""
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
-ninja          <true|false>    : set to true to build with ninja instead of make
-mpi            <true|false>    : set to true if you want to build the MPI version. Defaults to false.
-cmd            <true|false>    : set to true if you want to build RevStudio with GTK2+. Defaults to false.
-jupyter        <true|false>    : set to true if you want to build the jupyter version. Defaults to false.
-help           <true|false>    : update the help database and build the YAML help generator. Defaults to false.
-boost_root     string          : specify directory containing Boost headers (e.g. `/usr/include`). Defaults to unset.
-boost_lib      string          : specify directory containing Boost libraries. (e.g. `/usr/lib`). Defaults to unset.
-static_boost	<true|false>    : link using static Boost libraries. Defaults to false.
-j              integer         : the number of threads to use when compiling RevBayes. Defaults to 4.

You can also specify cmake variables as -DCMAKE_VAR1=value1 -DCMAKE_VAR2=value2

Examples:
  ./build.sh -mpi true -help true
  ./build.sh -boost_root /home/santa/boost_1.72
  ./build.sh -DBOOST_ROOT=/home/santa/boost_1.72
  ./build.sh -mpi true -DHELP=ON -DBOOST_ROOT=/home/santa/boost_1.72'
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

if [ "$mpi" = "true" ] ; then
    BUILD_DIR="build-mpi"
else
    BUILD_DIR="build"
fi

if [ -z "${exec_name}" ] ; then
    if [ "$mpi" = "true" ] ; then
        exec_name="rb-mpi"
    else
        exec_name=rb
    fi
fi

if [ "$travis" = "true" ]; then
    BUILD_DIR="build"
    CC=${C_COMPILER}
    CXX=${CXX_COMPILER}
    all_args="-travis true -mpi ${USE_MPI} -help true -exec_name rb"
    exec_name=rb
    help=true
fi

if [ "$debug" = "true" ] ; then
    cmake_args="-DCMAKE_BUILD_TYPE=DEBUG $cmake_args"
fi

if [ "$ninja" = "true" ] ; then
    cmake_args="$cmake_args -G Ninja"
fi

if [ "$mpi" = "true" ] ; then
    cmake_args="-DMPI=ON $cmake_args"
fi

if [ "$jupyter" = "true" ] ; then
    cmake_args="-DJUPYTER=ON $cmake_args"
fi

if [ "$cmd" = "true" ] ; then
    cmake_args="-DCMD_GTK=ON $cmake_args"
fi

if [ "$travis" = "true" ] ; then
    cmake_args="-DCONTINUOUS_INTEGRATION=TRUE $cmake_args"
fi

if [ -n "$boost_root" ] ; then
    cmake_args="-DBOOST_ROOT=\"${boost_root}\" $cmake_args"
fi

if [ -n "$boost_lib" ] ; then
    cmake_args="-DBOOST_LIBRARYDIR=\"${boost_lib}\" $cmake_args"
fi

if [ "$static_boost" = "true" ] ; then
    cmake_args="-DSTATIC_BOOST=ON $cmake_args"
fi

if [ "$help" = "true" ] ; then
    cmake_args="-DHELP=ON $cmake_args"
fi

echo "Building ${exec_name}"
cmake_args="-DRB_EXEC_NAME=${exec_name} $cmake_args"

# Just print the values of CC and CXX in case the user has some weird
# values set.
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

    ######### Generate help database
    if [ "$help" = "true" ]
    then
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
    if [ "$ninja" = "true" ] ; then
        echo "Running 'ninja -j $j' in $(pwd)"
        ninja -j $j
    else
        echo "Running 'make -j $j' in $(pwd)"
        make -j $j
    fi
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
