#!/bin/sh

set -e

#################
# command line options
# set default values
travis="false"
mpi="false"
help="false"

# parse command line arguments
while echo $1 | grep ^- > /dev/null; do
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

######## Generate CMakeLists.txt for subdirs
generate_subdir_cmake()
{
    local subdir="$1"
    local libname=$2
    if [ ! -d "$BUILD_DIR/$subdir" ]; then
        mkdir "$BUILD_DIR/$subdir"
    fi
    echo "set(${subdir}_FILES" > "$BUILD_DIR/$subdir/CMakeLists.txt"
    ### Treat all files in each subdir as source files
    find "$subdir" | grep -v "svn" | sed 's|^|${PROJECT_SOURCE_DIR}/|g' >> "$BUILD_DIR/$subdir/CMakeLists.txt"
    echo ")
add_library($libname \${${subdir}_FILES})"  >> "$BUILD_DIR/$subdir/CMakeLists.txt"
}

generate_subdir_cmake "libs" "rb-libs"

generate_subdir_cmake "core" "rb-core"

generate_subdir_cmake "revlanguage" "rb-parser"

# We will only USE this if we do subdir(help2yml)
generate_subdir_cmake "help2yml" "rb-help"

# We will only USE this if we do subdir(cmd)
generate_subdir_cmake "cmd" "rb-cmd-lib"
