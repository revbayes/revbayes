#!/bin/sh

set -e

SCRIPT_DIR=$(pwd)

if [ "$#" = 1 ] ; then
    BUILD_DIR="$1"
else
    BUILD_DIR=build
fi

if [ ! -e "$BUILD_DIR" ] ; then
   mkdir "$BUILD_DIR"
fi

BUILD_DIR=$(cd "$BUILD_DIR"; pwd)

SRC_DIR="$BUILD_DIR"/../../../src
cd "$BUILD_DIR"/../../../src

######### Generate generated_include_dirs.cmake ############
cat "$SCRIPT_DIR/cmake-fragments/CMakeLists.txt" > "$BUILD_DIR/CMakeLists.txt"

######### Generate generated_include_dirs.cmake ############
(
    echo '
    # TODO Split these up based on sub-package dependency
INCLUDE_DIRECTORIES('
    find libs core revlanguage -type d | grep -v "svn" | sed 's|^|    ${PROJECT_SOURCE_DIR}/|g' 
    echo ' ${Boost_INCLUDE_DIR} )
'
) > "$BUILD_DIR/generated_include_dirs.cmake"

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
