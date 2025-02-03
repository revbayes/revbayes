#!/bin/sh

# NOTE: All configuration is now done via cmake.
#
#       * Configuration options to this file are translated into cmake variables and passed to cmake.
#       * Options like -DKEY=VALUE are passed to cmake.
#
#       To debug configuration problems:
#         * look at src/CMakeLists.txt
#         * note what options this script passes to cmake (there's a log message)

# NOTE: Overview of what this script does:
# 1. Read command line flags
# 2. Generate cmake variables from command-line flags.
# 3. Create the build/ directory (if missing).
# 4. Update the version number            --> src/revlanguage/utils/GitVersion.cpp
# 5. Update the help database --> src/core/help/RbHelpDatabase.cpp
# 6. Run ./regenerate.sh
# 7. Run cmake <--- This is where the configuration actually happens
# 8. Run make or ninja to do the build.
# 9. Restore GitVersion.cpp from backup.

# If you change this script, please update the list above.


set -e

all_args="$@"

#################
# command line options
# set default values
debug="false"
mpi="false"
cmd="false"
boost_root=""
boost_lib=""
boost_include=""
boost_verbose=""
boost_debug=""
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
-boost_root     string          : specify directory containing Boost headers and libraries (e.g. `/usr/`). Defaults to unset.
-boost_lib      string          : specify directory containing Boost libraries. (e.g. `/usr/lib`). Defaults to unset.
-boost_include  string          : specify directory containing Boost libraries. (e.g. `/usr/include`). Defaults to unset.
-boost_verbose  <true|false>    : log some info about finding Boost
-boost_debug    <true|false>    : log MORE info about finding Boost
-static_boost	<true|false>    : link using static Boost libraries. Defaults to false.
-j              integer         : the number of threads to use when compiling RevBayes. Defaults to 4.

You can also specify cmake variables as -DCMAKE_VAR1=value1 -DCMAKE_VAR2=value2

Examples:
  ./build.sh -mpi true -help2yml true
  ./build.sh -boost_root /home/santa/installed-boost-1.72
  ./build.sh -boost_include /home/santa/boost_1_72_0/ -boost_lib /home/santa/boost_1_72_0/stage/lib
  ./build.sh -DBOOST_ROOT=/home/santa/installed-boost_1.72
  ./build.sh -mpi true -DBOOST_ROOT=/home/santa/installed-boost_1.72'
        exit
    fi

    case "$1" in -D*)
                     cmake_args="$cmake_args $1"
                     shift
                     continue
                     ;;
                 -*=*)
                     echo "$0: I don't understand '$1' - did you mean '$(echo $1 | sed 's/=/ /')'?"
                     exit 1
                     ;;
    esac

    # parse pairs
    eval $( echo $1 | sed 's/-//g' | tr -d '\012')=$2
    shift
    shift
done

if [ -z "${BUILD_DIR}" ] ; then
    if [ "$mpi" = "true" ] ; then
        BUILD_DIR="build-mpi"
    else
        BUILD_DIR="build"
    fi
fi

if [ -z "${exec_name}" ] ; then
    if [ "$mpi" = "true" ] ; then
        exec_name="rb-mpi"
    else
        exec_name=rb
    fi
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

if [ "$cmd" = "true" ] ; then
    cmake_args="-DCMD_GTK=ON $cmake_args"
fi

if [ -n "$jupyter" ] ; then
    echo "There is no longer a -jupyter <true|false> option to '$0'."
    echo "Jupyter functionality is now part of the standard rb application."
    echo
    echo "Run '$0 -h' to see available options."
    exit 1
fi

if [ -n "$boost_lib" ] && [ -n "$boost_include" ] ; then
    export BOOST_INCLUDEDIR="${boost_include}"
    export BOOST_LIBRARYDIR="${boost_lib}"
    unset BOOST_ROOT
    if [ -n "$boost_root" ] ; then
        echo "If you specify -boost_lib or -boost_include, then you cannot also specify -boost_root."
        exit 1
    fi
elif [ -n "$boost_lib" ] || [ -n "$boost_include" ] ; then
    echo "The flags -boost_lib and -boost_include must be given together"
    exit 1
elif [ -n "$boost_root" ] ; then
    export BOOST_ROOT="${boost_root}"
    unset BOOST_INCLUDEDIR
    unset BOOST_LIBRARYDIR
fi

if [ "$boost_verbose" = "true" ] ; then
    cmake_args="-DBoost_VERBOSE=ON $cmake_args"
fi

if [ "$boost_debug" = "true" ] ; then
    cmake_args="-DBoost_DEBUG=ON $cmake_args"
fi

if [ "$static_boost" = "true" ] ; then
    cmake_args="-DSTATIC_BOOST=ON $cmake_args"
fi

# generate rb-help2yml executable
# manually set DHELP=OFF to avoid
cmake_args="-DHELP=ON $cmake_args"

echo "RevBayes executable is '${exec_name}'"
cmake_args="-DRB_EXEC_NAME=${exec_name} $cmake_args"

if [ "$1" = "clean" ]
then
    rm -rf ${BUILD_DIR}
    exit 1
fi

if [ ! -d ${BUILD_DIR} ]; then
    mkdir ${BUILD_DIR}
fi

######## generate git version number
./generate_version_number.sh
echo " Saving old GitVersion.cpp in projects/cmake/GitVersion_backup.cpp"
if [ -e ../../src/revlanguage/utils/GitVersion.cpp ] ; then
    cp ../../src/revlanguage/utils/GitVersion.cpp GitVersion_backup.cpp
fi
echo " Copying current GitVersion.cpp to src/revlanguage/utils"
mv GitVersion.cpp ../../src/revlanguage/utils/


######### Generate help database
(
cd ../
./generate_help.sh
)

######## Generate some files for cmake
echo "Running './regenerate.sh $(pwd)/$BUILD_DIR"
./regenerate.sh $(pwd)/$BUILD_DIR
cd ${BUILD_DIR}
echo
echo


######### Print some environment variables
# * This can alert the user if some weird values have been set.
# * This also helps us replicate the call to cmake.
echo "Note these environment variables:"
for var in CC CXX CFLAGS CPPFLAGS CXXFLAGS LDFLAGS BOOST_ROOT BOOST_INCLUDEDIR BOOST_LIBRARYDIR ; do
    cmd="if [ -n \"\${$var}\" ] ; then echo \"  ${var}=\${$var}\"; fi"
    eval $cmd
done
echo

######### Actually run cmake
echo "Running 'cmake ../../../src $cmake_args' in $(pwd)"
cmake ../../../src $cmake_args
echo


######### Do the build
if [ "$ninja" = "true" ] ; then
    echo "Running 'ninja -j $j' in $(pwd)"
    ninja -j $j
else
    echo "Running 'make -j $j' in $(pwd)"
    make -j $j
fi
cd ..


####### Restore GitVersion.cpp from backup
if [ -e  GitVersion_backup.cpp ] ; then
    cp GitVersion_backup.cpp ../../src/revlanguage/utils/GitVersion.cpp
    rm GitVersion_backup.cpp
fi

