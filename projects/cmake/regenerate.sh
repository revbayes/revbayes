#!/bin/sh

set -e

#################
# command line options
# set default values
debug="false"
mac="false"
travis="false"
win="false"
mpi="false"
gentoo="false"
help="false"
jupyter="false"
boost_root=""
boost_lib=""

# parse command line arguments
while echo $1 | grep ^- > /dev/null; do
    # intercept help while parsing "-key value" pairs
    if [ "$1" = "--help" ] || [ "$1" = "-h" ]
    then
        echo '
The minimum steps to build RevBayes after running this script is:
cmake .
make

Command line options are:
-h                              : print this help and exit.
-mac            <true|false>    : set to true if you are building for a OS X - compatible with 10.6 and higher. Defaults to false.
-win            <true|false>    : set to true if you are building on a Windows system. Defaults to false.
-mpi            <true|false>    : set to true if you want to build the MPI version. Defaults to false.
-help           <true|false>    : Update the help database and build the YAML help generator. Defaults to false.
'
# secret test option
# -jupyter        <true|false>    : set to true if you want ot buikd the jupyter version. Defaults to false.
        exit
    fi

    # skip cmake args
    case "$1" in -D*)
                     shift
                     continue
                     ;;
    esac

    # parse pairs
    eval $( echo $1 | sed 's/-//g' | tr -d '\012')=$2
    shift
    shift
done



#################
# generate cmake configuration

if [ "${mpi}" = "true" ] && [ "${travis}" = "false" ]; then
    BUILD_DIR="$(pwd)/build-mpi"
    mkdir -p "${BUILD_DIR}"
    echo $BUILD_DIR
else
    BUILD_DIR="$(pwd)/build"
    mkdir -p "${BUILD_DIR}"
    echo $BUILD_DIR
fi

SCRIPT_DIR=$(pwd)

cd "$BUILD_DIR"/../../../src

######### Generate CMakeLists.txt: top #########################
cat "$SCRIPT_DIR/cmake-fragments/CMakeLists-top.txt" > "$BUILD_DIR/CMakeLists.txt"

######### Generate CMakeLists.txt: list of include directories #
echo '
# TODO Split these up based on sub-package dependency
INCLUDE_DIRECTORIES(' >> "$BUILD_DIR/CMakeLists.txt"
find libs core revlanguage -type d | grep -v "svn" | sed 's|^|    ${PROJECT_SOURCE_DIR}/|g' >> "$BUILD_DIR/CMakeLists.txt"
echo ' ${Boost_INCLUDE_DIR} )
'>> "$BUILD_DIR/CMakeLists.txt"

######### Generate CMakeLists.txt: bottom ######################
cat "$SCRIPT_DIR/cmake-fragments/CMakeLists-bottom.txt" >> "$BUILD_DIR/CMakeLists.txt"

######### Generate help database
if [ "$help" = "true" ]
then
    echo "Generating help database"
    perl ../help/md2help.pl ../help/md/*.md > core/help/RbHelpDatabase.cpp
fi

######### Generate libs/CMakeLists.txt #########################
if [ ! -d "$BUILD_DIR/libs" ]; then
    mkdir "$BUILD_DIR/libs"
fi
echo 'set(LIBS_FILES' > "$BUILD_DIR/libs/CMakeLists.txt"
find libs | grep -v "svn" | sed 's|^|${PROJECT_SOURCE_DIR}/|g' >> "$BUILD_DIR/libs/CMakeLists.txt"
echo ')
add_library(rb-libs ${LIBS_FILES})'  >> "$BUILD_DIR/libs/CMakeLists.txt"

######### Generate core/CMakeLists.txt #########################
if [ ! -d "$BUILD_DIR/core" ]; then
    mkdir "$BUILD_DIR/core"
fi
echo 'set(CORE_FILES' > "$BUILD_DIR/core/CMakeLists.txt"
find core | grep -v "svn" | sed 's|^|${PROJECT_SOURCE_DIR}/|g' >> "$BUILD_DIR/core/CMakeLists.txt"
echo ')
add_library(rb-core ${CORE_FILES})'  >> "$BUILD_DIR/core/CMakeLists.txt"

######### Generate revlanguage/CMakeLists.txt ##################
if [ ! -d "$BUILD_DIR/revlanguage" ]; then
    mkdir "$BUILD_DIR/revlanguage"
fi
echo 'set(PARSER_FILES' > "$BUILD_DIR/revlanguage/CMakeLists.txt"
find revlanguage | grep -v "svn" | sed 's|^|${PROJECT_SOURCE_DIR}/|g' >> "$BUILD_DIR/revlanguage/CMakeLists.txt"
echo ')
add_library(rb-parser ${PARSER_FILES})'  >> "$BUILD_DIR/revlanguage/CMakeLists.txt"

######### Generate help2yml/CMakeLists.txt #####################
# We will only USE this if we do subdir(help2yml)
if [ ! -d "$BUILD_DIR/help2yml" ]; then
    mkdir "$BUILD_DIR/help2yml"
fi
echo 'set(HELP_FILES' > "$BUILD_DIR/help2yml/CMakeLists.txt"
find help2yml | grep -v "svn" | sed 's|^|${PROJECT_SOURCE_DIR}/|g' >> "$BUILD_DIR/help2yml/CMakeLists.txt"
echo ')
add_library(rb-help ${HELP_FILES})'  >> "$BUILD_DIR/help2yml/CMakeLists.txt"

######### Generate cmd/CMakeLists.txt ##########################
# We will only USE this if we do subdir(cmd)
if [ ! -d "$BUILD_DIR/cmd" ]; then
    mkdir "$BUILD_DIR/cmd"
fi
echo 'SET(CMD_FILES' > "$BUILD_DIR/cmd/CMakeLists.txt"
find cmd | grep -v "svn" | sed 's|^|${PROJECT_SOURCE_DIR}/|g' >> "$BUILD_DIR/cmd/CMakeLists.txt"
echo ')
ADD_LIBRARY(rb-cmd-lib ${CMD_FILES})'  >> "$BUILD_DIR/cmd/CMakeLists.txt"



