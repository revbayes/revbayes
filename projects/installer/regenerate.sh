#!/bin/sh

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
boost_root="${MY_BOOST_ROOT}"
boost_lib="${MY_BOOST_ROOT}/lib"
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

cd "$BUILD_DIR"/../../../src

echo 'cmake_minimum_required(VERSION 2.6)
project(RevBayes)

# Consider:
# -Wno-sign-compare
# -D_GLIBCXX_DEBUG

# This is the RIGHT way, but requires cmake version >=3:
#   set(CMAKE_CXX_STANDARD 11)
# RHEL 7 compute clusters may have cmake 2.8.12
#
# So, we add the flag directly instead.
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

if(NOT (${CMAKE_VERSION} VERSION_LESS "2.8.0"))
  find_program(CCACHE_PROGRAM ccache)
  if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
  endif()
endif()

' > "$BUILD_DIR/CMakeLists.txt"


if [ "${exec_name}" = "" ]; then
    if [ "${mpi}" = "true" ]; then
        exec_name="rb-mpi"
    else
        exec_name="rb"
    fi
fi

echo "Building ${exec_name}"

echo "set (RB_EXEC_NAME \"${exec_name}\")" >> "$BUILD_DIR/CMakeLists.txt"
echo "set (LOCAL_BOOST_ROOT \"${boost_root}\")" >> "$BUILD_DIR/CMakeLists.txt"
echo "set (LOCAL_BOOST_LIBRARY \"${boost_lib}\")" >> "$BUILD_DIR/CMakeLists.txt"

if [ "$debug" = "true" ]
then
echo '
# -Woverloaded-virtual has some false-positives with GCC
# We should ultimiately remove -Wno-reorder -Wno-unused-variable -Wno-unused-but-set-variable
# But there are so many of them we cant see the really bad warnings.
# So, disable those warnings for now.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -Wall -msse -msse2 -msse3")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O0 -Wall")
'  >> "$BUILD_DIR/CMakeLists.txt"
elif [ "$mac" = "true" ]
then
echo '
set(CMAKE_OSX_DEPLOYMENT_TARGET "10.6")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -msse -msse2 -msse3")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
'  >> "$BUILD_DIR/CMakeLists.txt"
elif [ "$win" = "true" ]
then
echo '
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -msse -msse2 -msse3")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
'  >> "$BUILD_DIR/CMakeLists.txt"
else
echo '
if (CMAKE_SYSTEM_PROCESSOR MATCHES "^arm*|aarch64")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
   add_definitions(-DRB_ARM)
else()
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -msse -msse2 -msse3")
endif()
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
'  >> "$BUILD_DIR/CMakeLists.txt"
fi

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

echo '
# Add extra CMake libraries into ./CMake
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake ${CMAKE_MODULE_PATH})

# Set source root relate to project file
set(PROJECT_SOURCE_DIR ${CMAKE_SOURCE_DIR}/../../../src)

MESSAGE("My Boost information:")
MESSAGE("  Boost_INCLUDE_DIRS: ${LOCAL_BOOST_ROOT}")
MESSAGE("  Boost_LIBRARIES: ${LOCAL_BOOST_LIBRARY}")

if ( NOT ${LOCAL_BOOST_ROOT} STREQUAL "")

   SET(BOOST_ROOT ${LOCAL_BOOST_ROOT})
   SET(BOOST_INCLUDEDIR ${LOCAL_BOOST_ROOT}/include)
   SET(BOOST_LIBRARYDIR ${LOCAL_BOOST_ROOT}/lib)
   SET(Boost_NO_SYSTEM_PATHS ON)
   SET(Boost_USE_STATIC_LIBS ON)

endif()

find_package(Boost
1.60.0
COMPONENTS regex
program_options
thread
system
filesystem
date_time
serialization REQUIRED)
MESSAGE("Boost information:")
MESSAGE("  Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIR}")
MESSAGE("  Boost_LIBRARIES: ${Boost_LIBRARIES}")
MESSAGE("  Boost_LIBRARY_DIRS: ${Boost_LIBRARY_DIRS}")
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})


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
perl ../help/md2help.pl ../help/md/* > core/help/RbHelpDatabase.cpp

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
echo '

# Use the package PkgConfig to detect GTK+ headers/library files
FIND_PACKAGE(PkgConfig REQUIRED)
PKG_CHECK_MODULES(GTK REQUIRED gtk+-2.0)
#PKG_CHECK_MODULES(GTK REQUIRED gtk+-3.0)
' >> $BUILD_DIR/CMakeLists.txt

# Setup CMake to use GTK+, tell the compiler where to look for headers
# and to the linker where to look for libraries
echo '
INCLUDE_DIRECTORIES(${GTK_INCLUDE_DIRS})
LINK_DIRECTORIES(${GTK_LIBRARY_DIRS})
' >> $BUILD_DIR/CMakeLists.txt

echo '
# Add other flags to the compiler
ADD_DEFINITIONS(${GTK_CFLAGS_OTHER})

# Add an executable compiled from hello.c
ADD_EXECUTABLE(RevStudio ${PROJECT_SOURCE_DIR}/cmd/main.cpp)
' >> $BUILD_DIR/CMakeLists.txt

# Link the target to the GTK+ libraries
echo '
TARGET_LINK_LIBRARIES(RevStudio rb-cmd-lib rb-parser rb-core rb-libs ${Boost_LIBRARIES} ${GTK_LIBRARIES})
' >> $BUILD_DIR/CMakeLists.txt

echo '
SET_TARGET_PROPERTIES(RevStudio PROPERTIES PREFIX "../")

add_subdirectory(cmd)
' >> $BUILD_DIR/CMakeLists.txt

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

target_link_libraries(${RB_EXEC_NAME} rb-parser rb-core rb-libs ${Boost_LIBRARIES} ${CMAKE_DL_LIBS})

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
