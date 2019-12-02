#!/bin/sh
HERE="$(pwd)/build"
mkdir -p "${HERE}"
echo $HERE

#################
# command line options
# set default values
boost="true"
debug="false"
mac="false"
travis="false"
win="false"
mpi="false"
gentoo="false"
help="false"
jupyter="false"

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
-boost          <true|false>    : true (re)compiles boost libs, false dont. Defaults to true.
-mac            <true|false>    : set to true if you are building for a OS X - compatible with 10.6 and higher. Defaults to false.
-win            <true|false>    : set to true if you are building on a Windows system. Defaults to false.
-mpi            <true|false>    : set to true if you want to build the MPI version. Defaults to false.
-help           <true|false>    : Build the help generator. Defaults to false.
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
# build boost libraries separately

if [ "$boost" = "true" ]
then
    echo 'Building boost libraries'
    echo 'you can turn this off with argument "-boost false"'

    cd ../../boost_1_60_0
    rm ./project-config.jam*  # clean up from previous runs

    if [ "$mac" = "true" ]
    then
        ./bootstrap.sh --with-libraries=atomic,chrono,filesystem,system,regex,thread,date_time,program_options,math,serialization
        ./b2 link=static
    elif [ "$win" = "true" ]
    then
        ./bootstrap.sh --with-libraries=atomic,chrono,filesystem,system,regex,thread,date_time,program_options,math,serialization --with-toolset=mingw
        ./b2 link=static
    elif [ "$gentoo" = "true" ]
    then
        ./bootstrap.sh --with-libraries=atomic,chrono,filesystem,system,regex,thread,date_time,program_options,math,serialization
        ./b2 link=static --ignore-site-config
    else
        ./bootstrap.sh --with-libraries=atomic,chrono,filesystem,system,regex,thread,date_time,program_options,math,serialization
        ./b2 link=static
    fi

else
    echo 'not building boost libraries'
fi


#################
# generate cmake configuration
cd "$HERE"/../
cd ../../src

echo 'cmake_minimum_required(VERSION 2.6)
project(RevBayes)

# Consider:
# -Wno-sign-compare
# -D_GLIBCXX_DEBUG

# Default compiler flags
#if (WIN32)
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -g -pg -static -msse -msse2 -msse3 -Wall -Wno-sign-compare")
#    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wall -g -pg -static")
#else ()
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -g -pg -msse -msse2 -msse3 -stdlib=libstdc++ -Wall -Wno-sign-compare")
#    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wall -g -pg")
#endif ()

' > "$HERE/CMakeLists.txt"

if [ "$debug" = "true" ]
then
echo '
# -Woverloaded-virtual has some false-positives with GCC
# We should ultimiately remove -Wno-reorder -Wno-unused-variable -Wno-unused-but-set-variable
# But there are so many of them we cant see the really bad warnings.
# So, disable those warnings for now.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -Wall -msse -msse2 -msse3")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O0 -Wall")
'  >> "$HERE/CMakeLists.txt"
elif [ "$mac" = "true" ]
then
echo '
set(CMAKE_OSX_DEPLOYMENT_TARGET "10.6")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -msse -msse2 -msse3")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
'  >> "$HERE/CMakeLists.txt"
elif [ "$win" = "true" ]
then
echo '
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -msse -msse2 -msse3 -static -std=gnu++98")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -static")
'  >> "$HERE/CMakeLists.txt"
else
echo '
if (CMAKE_SYSTEM_PROCESSOR MATCHES "^arm*|aarch64")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
   add_definitions(-DRB_ARM)
else()
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -msse -msse2 -msse3")
endif()
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
'  >> "$HERE/CMakeLists.txt"
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
'  >> "$HERE/CMakeLists.txt"
fi

if [ "$win" = "true" ]
then
echo '
add_definitions(-DRB_WIN)
'  >> "$HERE/CMakeLists.txt"
fi


if [ "$jupyter" = "true" ]
then
echo "JUPYTER!"
echo '
add_definitions(-DRB_XCODE)
'  >> "$HERE/CMakeLists.txt"
fi

if [ "$travis" = "true" ]
then
    echo 'set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g0 -O2")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g0 -O2")
' >> "$HERE/CMakeLists.txt"
fi

echo '
# Add extra CMake libraries into ./CMake
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake ${CMAKE_MODULE_PATH})

# Set source root relate to project file
set(PROJECT_SOURCE_DIR ${CMAKE_SOURCE_DIR}/../../../src)

option(INTERNAL_BOOST "Use the version of boost shipped with revbayes" ON)

if (INTERNAL_BOOST)
   SET(BOOST_ROOT "../../../boost_1_60_0")
   SET(BOOST_LIBRARY "../../../boost_1_60_0/stage/lib")
   SET(Boost_NO_SYSTEM_PATHS ON)
   SET(Boost_USE_STATIC_RUNTIME ON)
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
INCLUDE_DIRECTORIES(' >> "$HERE/CMakeLists.txt"
find libs core revlanguage help -type d | grep -v "svn" | sed 's|^|    ${PROJECT_SOURCE_DIR}/|g' >> "$HERE/CMakeLists.txt"
echo ' ${Boost_INCLUDE_DIR} )


####

# Split into much smaller libraries
add_subdirectory(libs)
add_subdirectory(core)
add_subdirectory(revlanguage)


############# executables #################
# basic rev-bayes binary
' >> $HERE/CMakeLists.txt

if [ "$help" = "true" ]
then
echo '
add_subdirectory(help)
' >> $HERE/CMakeLists.txt
fi

if [ "$help" = "true" ]
then
echo '
add_executable(rb-help ${PROJECT_SOURCE_DIR}/help/YAMLHelpGenerator.cpp)

target_link_libraries(rb-help rb-parser rb-core libs help_yaml ${Boost_LIBRARIES})
set_target_properties(rb-help PROPERTIES PREFIX "../")
' >> $HERE/CMakeLists.txt
elif [ "$jupyter" = "true" ]
then
echo "more jupyter!"
echo '
add_executable(rb-jupyter ${PROJECT_SOURCE_DIR}/revlanguage/main.cpp)

target_link_libraries(rb-jupyter rb-parser rb-core libs ${Boost_LIBRARIES})
set_target_properties(rb-jupyter PROPERTIES PREFIX "../")
' >> $HERE/CMakeLists.txt
elif [ "$cmd" = "true" ]
then
echo '

# Use the package PkgConfig to detect GTK+ headers/library files
FIND_PACKAGE(PkgConfig REQUIRED)
PKG_CHECK_MODULES(GTK REQUIRED gtk+-2.0)
#PKG_CHECK_MODULES(GTK REQUIRED gtk+-3.0)
' >> $HERE/CMakeLists.txt

# Setup CMake to use GTK+, tell the compiler where to look for headers
# and to the linker where to look for libraries
echo '
INCLUDE_DIRECTORIES(${GTK_INCLUDE_DIRS})
LINK_DIRECTORIES(${GTK_LIBRARY_DIRS})
' >> $HERE/CMakeLists.txt

echo '
# Add other flags to the compiler
ADD_DEFINITIONS(${GTK_CFLAGS_OTHER})

# Add an executable compiled from hello.c
ADD_EXECUTABLE(RevStudio ${PROJECT_SOURCE_DIR}/cmd/main.cpp)
' >> $HERE/CMakeLists.txt

# Link the target to the GTK+ libraries
echo '
TARGET_LINK_LIBRARIES(RevStudio rb-cmd-lib rb-parser rb-core libs ${Boost_LIBRARIES} ${GTK_LIBRARIES})
' >> $HERE/CMakeLists.txt

echo '
SET_TARGET_PROPERTIES(RevStudio PROPERTIES PREFIX "../")

add_subdirectory(cmd)
' >> $HERE/CMakeLists.txt

if [ ! -d "$HERE/cmd" ]; then
mkdir "$HERE/cmd"
fi
echo 'SET(CMD_FILES' > "$HERE/cmd/CMakeLists.txt"
find cmd | grep -v "svn" | sed 's|^|${PROJECT_SOURCE_DIR}/|g' >> "$HERE/cmd/CMakeLists.txt"
echo ')
ADD_LIBRARY(rb-cmd-lib ${CMD_FILES})'  >> "$HERE/cmd/CMakeLists.txt"

else
echo '
add_executable(rb ${PROJECT_SOURCE_DIR}/revlanguage/main.cpp)

target_link_libraries(rb rb-parser rb-core libs ${Boost_LIBRARIES})

set_target_properties(rb PROPERTIES PREFIX "../")
' >> $HERE/CMakeLists.txt
if [ "$mpi" = "true" ] ; then
    echo 'target_link_libraries(rb ${MPI_LIBRARIES})
' >> $HERE/CMakeLists.txt
fi
fi

if [ ! -d "$HERE/libs" ]; then
mkdir "$HERE/libs"
fi
echo 'set(LIBS_FILES' > "$HERE/libs/CMakeLists.txt"
find libs | grep -v "svn" | sed 's|^|${PROJECT_SOURCE_DIR}/|g' >> "$HERE/libs/CMakeLists.txt"
echo ')
add_library(libs ${LIBS_FILES})'  >> "$HERE/libs/CMakeLists.txt"

if [ ! -d "$HERE/help" ]; then
mkdir "$HERE/help"
fi
echo 'set(HELP_FILES' > "$HERE/help/CMakeLists.txt"
find help | grep -v "svn" | sed 's|^|${PROJECT_SOURCE_DIR}/|g' >> "$HERE/help/CMakeLists.txt"
echo ')
add_library(help_yaml ${HELP_FILES})'  >> "$HERE/help/CMakeLists.txt"

if [ ! -d "$HERE/core" ]; then
mkdir "$HERE/core"
fi
echo 'set(CORE_FILES' > "$HERE/core/CMakeLists.txt"
find core | grep -v "svn" | sed 's|^|${PROJECT_SOURCE_DIR}/|g' >> "$HERE/core/CMakeLists.txt"
echo ')
add_library(rb-core ${CORE_FILES})'  >> "$HERE/core/CMakeLists.txt"

if [ ! -d "$HERE/revlanguage" ]; then
mkdir "$HERE/revlanguage"
fi
echo 'set(PARSER_FILES' > "$HERE/revlanguage/CMakeLists.txt"
find revlanguage | grep -v "svn" | sed 's|^|${PROJECT_SOURCE_DIR}/|g' >> "$HERE/revlanguage/CMakeLists.txt"
echo ')
add_library(rb-parser ${PARSER_FILES})'  >> "$HERE/revlanguage/CMakeLists.txt"
