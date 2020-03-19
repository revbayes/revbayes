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
exec_name=""

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

cat "$SCRIPT_DIR/cmake-fragments/CMakeLists-top.txt" > "$BUILD_DIR/CMakeLists.txt"

if [ "${exec_name}" = "" ]; then
    if [ "${mpi}" = "true" ]; then
        exec_name="rb-mpi"
    else
        exec_name="rb"
    fi
fi

echo "set (LOCAL_BOOST_ROOT \"${boost_root}\")" >> "$BUILD_DIR/CMakeLists.txt"
echo "set (LOCAL_BOOST_LIBRARY \"${boost_lib}\")" >> "$BUILD_DIR/CMakeLists.txt"

if [ "$mpi" = "true" ]
then
    echo '
add_definitions(-DRB_MPI)
#add_definitions(-DDEBUG_MPI_MCA)
# Require MPI for this project:
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS}")
'  >> "$BUILD_DIR/CMakeLists.txt"
fi

if [ "$win" = "true" ]
then
    echo '
add_definitions(-DRB_WIN)
'  >> "$BUILD_DIR/CMakeLists.txt"
fi


if [ "$jupyter" = "true" ]
then
echo "JUPYTER!"
echo '
add_definitions(-DRB_XCODE)
'  >> "$BUILD_DIR/CMakeLists.txt"
fi

if [ "$travis" = "true" ]
then
    echo 'set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g0 -O2")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g0 -O2")
' >> "$BUILD_DIR/CMakeLists.txt"
fi

## Emit CMakeLists-boost.txt
cat "$SCRIPT_DIR/cmake-fragments/CMakeLists-boost.txt" >> "$BUILD_DIR/CMakeLists.txt"

echo '
# TODO Split these up based on sub-package dependency
INCLUDE_DIRECTORIES(' >> "$BUILD_DIR/CMakeLists.txt"
find libs core revlanguage -type d | grep -v "svn" | sed 's|^|    ${PROJECT_SOURCE_DIR}/|g' >> "$BUILD_DIR/CMakeLists.txt"
echo ' ${Boost_INCLUDE_DIR} )


####

# Split into much smaller libraries
add_subdirectory(libs)
add_subdirectory(core)
add_subdirectory(revlanguage)


############# executables #################
# basic rev-bayes binary
' >> $BUILD_DIR/CMakeLists.txt

if [ "$help" = "true" ]
then

    #################
    # generate help database
    echo "Generating help database"
    perl ../help/md2help.pl ../help/md/*.md > core/help/RbHelpDatabase.cpp

    echo '
add_subdirectory(help2yml)
' >> $BUILD_DIR/CMakeLists.txt

    echo '
add_executable(${RB_EXEC_NAME}-help2yml ${PROJECT_SOURCE_DIR}/help2yml/main.cpp)

target_link_libraries(${RB_EXEC_NAME}-help2yml rb-help rb-parser rb-core rb-libs rb-parser ${Boost_LIBRARIES})
set_target_properties(${RB_EXEC_NAME}-help2yml PROPERTIES PREFIX "../")
' >> $BUILD_DIR/CMakeLists.txt

    if [ ! -d "$BUILD_DIR/help2yml" ]; then
        mkdir "$BUILD_DIR/help2yml"
    fi
    echo 'set(HELP_FILES' > "$BUILD_DIR/help2yml/CMakeLists.txt"
    find help2yml | grep -v "svn" | sed 's|^|${PROJECT_SOURCE_DIR}/|g' >> "$BUILD_DIR/help2yml/CMakeLists.txt"
    echo ')
add_library(rb-help ${HELP_FILES})'  >> "$BUILD_DIR/help2yml/CMakeLists.txt"

    if [ "$mpi" = "true" ] ; then
        echo 'target_link_libraries(${RB_EXEC_NAME}-help2yml ${MPI_LIBRARIES})
' >> $BUILD_DIR/CMakeLists.txt
    fi
fi

if [ "$jupyter" = "true" ]
then
    echo "more jupyter!"
    echo '
add_executable(rb-jupyter ${PROJECT_SOURCE_DIR}/revlanguage/main.cpp)

target_link_libraries(rb-jupyter rb-parser rb-core rb-libs ${Boost_LIBRARIES})
set_target_properties(rb-jupyter PROPERTIES PREFIX "../")
' >> $BUILD_DIR/CMakeLists.txt
elif [ "$cmd" = "true" ]
then
    cat "$SCRIPT_DIR/cmake-fragments/CMakeLists-RevStudio.txt" >> "$BUILD_DIR/CMakeLists.txt"

    if [ ! -d "$BUILD_DIR/cmd" ]; then
        mkdir "$BUILD_DIR/cmd"
    fi
    echo 'SET(CMD_FILES' > "$BUILD_DIR/cmd/CMakeLists.txt"
    find cmd | grep -v "svn" | sed 's|^|${PROJECT_SOURCE_DIR}/|g' >> "$BUILD_DIR/cmd/CMakeLists.txt"
    echo ')
ADD_LIBRARY(rb-cmd-lib ${CMD_FILES})'  >> "$BUILD_DIR/cmd/CMakeLists.txt"

else
    echo '
add_executable(${RB_EXEC_NAME} ${PROJECT_SOURCE_DIR}/revlanguage/main.cpp)

target_link_libraries(${RB_EXEC_NAME} rb-parser rb-core rb-libs ${Boost_LIBRARIES})

set_target_properties(${RB_EXEC_NAME} PROPERTIES PREFIX "../")
' >> $BUILD_DIR/CMakeLists.txt
    if [ "$mpi" = "true" ] ; then
        echo 'target_link_libraries(${RB_EXEC_NAME} ${MPI_LIBRARIES})
' >> $BUILD_DIR/CMakeLists.txt
    fi
fi

if [ ! -d "$BUILD_DIR/libs" ]; then
    mkdir "$BUILD_DIR/libs"
fi
echo 'set(LIBS_FILES' > "$BUILD_DIR/libs/CMakeLists.txt"
find libs | grep -v "svn" | sed 's|^|${PROJECT_SOURCE_DIR}/|g' >> "$BUILD_DIR/libs/CMakeLists.txt"
echo ')
add_library(rb-libs ${LIBS_FILES})'  >> "$BUILD_DIR/libs/CMakeLists.txt"

if [ ! -d "$BUILD_DIR/core" ]; then
    mkdir "$BUILD_DIR/core"
fi
echo 'set(CORE_FILES' > "$BUILD_DIR/core/CMakeLists.txt"
find core | grep -v "svn" | sed 's|^|${PROJECT_SOURCE_DIR}/|g' >> "$BUILD_DIR/core/CMakeLists.txt"
echo ')
add_library(rb-core ${CORE_FILES})'  >> "$BUILD_DIR/core/CMakeLists.txt"

if [ ! -d "$BUILD_DIR/revlanguage" ]; then
    mkdir "$BUILD_DIR/revlanguage"
fi
echo 'set(PARSER_FILES' > "$BUILD_DIR/revlanguage/CMakeLists.txt"
find revlanguage | grep -v "svn" | sed 's|^|${PROJECT_SOURCE_DIR}/|g' >> "$BUILD_DIR/revlanguage/CMakeLists.txt"
echo ')
add_library(rb-parser ${PARSER_FILES})'  >> "$BUILD_DIR/revlanguage/CMakeLists.txt"

cat "$SCRIPT_DIR/cmake-fragments/CMakeLists-bottom.txt" >> "$BUILD_DIR/CMakeLists.txt"
