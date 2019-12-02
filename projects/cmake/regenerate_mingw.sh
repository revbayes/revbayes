#!/bin/sh
HERE=$(pwd)/build
echo $HERE

#################
# command line options
# set default values
boost="true"
debug="false"
mac="false"
win="false"
mpi="false"
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
-h | -help                      : print this help and exit.
-boost          <true|false>    : true (re)compiles boost libs, false dont. Defaults to true.
-mac            <true|false>    : set to true if you are building for a OS X - compatible with 10.6 and higher. Defaults to false.
-win            <true|false>    : set to true if you are building on a Windows system. Defaults to false.
-mpi            <true|false>    : set to true if you want to build the MPI version. Defaults to false.
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
./bootstrap.sh --with-libraries=atomic,chrono,filesystem,system,regex,thread,date_time,program_options,math,serialization,signals
./b2 link=static
else
./bootstrap.sh --with-libraries=atomic,chrono,filesystem,system,regex,thread,date_time,program_options,math,serialization,signals --with-toolset=mingw
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
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -g -pg -static -msse -msse2 -msse3")
#    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wall -g -pg -static")
#else ()
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -g -pg -msse -msse2 -msse3 -stdlib=libstdc++")
#    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wall -g -pg")
#endif ()

' > "$HERE/CMakeLists.txt"

if [ "$debug" = "true" ]
then
echo '
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
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -msse -msse2 -msse3")
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
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
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

echo '
# Add extra CMake libraries into ./CMake
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake ${CMAKE_MODULE_PATH})

# Set source root relate to project file
set(PROJECT_SOURCE_DIR ${CMAKE_SOURCE_DIR}/../../../src)



SET(BOOST_ROOT ../../../boost_1_60_0)
SET(BOOST_LIBRARY ../../../boost_1_60_0/stage/lib)
SET(Boost_USE_STATIC_RUNTIME true)
SET(Boost_USE_STATIC_LIBS ON)
find_package(Boost
1.60.0
COMPONENTS regex
program_options
thread
signals
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
find libs core revlanguage -type d | grep -v "svn" | sed 's|^|    ${PROJECT_SOURCE_DIR}/|g' >> "$HERE/CMakeLists.txt"
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
add_subdirectory(revlanguage)
' >> $HERE/CMakeLists.txt
fi

if [ "$mpi" = "true" ]
then
echo '
add_executable(rb-mpi ${PROJECT_SOURCE_DIR}/revlanguage/main.cpp)

target_link_libraries(rb-mpi rb-parser rb-core libs ${Boost_LIBRARIES} ${MPI_LIBRARIES})
set_target_properties(rb-mpi PROPERTIES PREFIX "../")
' >> $HERE/CMakeLists.txt
elif [ "$help" = "true" ]
then
echo '
add_executable(rb ${PROJECT_SOURCE_DIR}/tool/help/HtmlHelpGenerator.cpp)

target_link_libraries(rb rb-parser rb-core libs ${Boost_LIBRARIES})
set_target_properties(rb PROPERTIES PREFIX "../")
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

# Setup CMake to use GTK+, tell the compiler where to look for headers
# and to the linker where to look for libraries
INCLUDE_DIRECTORIES( /mingw64/include/gtk-2.0;/mingw64/lib/gtk-2.0/include;/mingw64/include/pango-1.0;/mingw64/include/fribidi;/mingw64/include/cairo;/mingw64/include/atk-1.0;/mingw64/include/cairo;/mingw64/include/pixman-1;/mingw64/include;/mingw64/include/freetype2;/mingw64/include;/mingw64/include/harfbuzz;/mingw64/include/libpng16;/mingw64/include/gdk-pixbuf-2.0;/mingw64/include/libpng16;/mingw64/include;/mingw64/include/glib-2.0;/mingw64/lib/glib-2.0/include;/mingw64/include )
LINK_DIRECTORIES( /mingw64/lib )

# Add other flags to the compiler
ADD_DEFINITIONS(${GTK_CFLAGS_OTHER})

# Add an executable compiled from hello.c
ADD_EXECUTABLE(RevStudio ${PROJECT_SOURCE_DIR}/cmd/main.cpp)

# Link the target to the GTK+ libraries
TARGET_LINK_LIBRARIES(RevStudio rb-cmd-lib rb-parser rb-core libs ${Boost_LIBRARIES} 
"/mingw64/lib/libgtk-win32-2.0.dll.a" 
"/mingw64/lib/libgdk-win32-2.0.dll.a" 
"/mingw64/lib/libpangowin32-1.0.dll.a" 
"/mingw64/lib/libpangocairo-1.0.dll.a" 
"/mingw64/lib/libpango-1.0.dll.a" 
"/mingw64/lib/libfribidi.dll.a" 
"/mingw64/lib/libatk-1.0.dll.a" 
"/mingw64/lib/libcairo.dll.a" 
"/mingw64/lib/libgdk_pixbuf-2.0.dll.a" 
"/mingw64/lib/libgio-2.0.dll.a" 
"/mingw64/lib/libgobject-2.0.dll.a" 
"/mingw64/lib/libglib-2.0.dll.a" 
"/mingw64/lib/libintl.dll.a"
)
SET_TARGET_PROPERTIES(RevStudio PROPERTIES PREFIX "../")
#SET_PROPERTY(TARGET rb-cmd PROPERTY CXX_STANDARD 98)

add_subdirectory(cmd)
' >> $HERE/CMakeLists.txt

if [ ! -d "$HERE/cmd" ]; then
mkdir "$HERE/cmd"
fi
echo 'set(CMD_FILES' > "$HERE/cmd/CMakeLists.txt"
find cmd | grep -v "svn" | sed 's|^|${PROJECT_SOURCE_DIR}/|g' >> "$HERE/cmd/CMakeLists.txt"
echo ')
ADD_LIBRARY(rb-cmd-lib ${CMD_FILES})'  >> "$HERE/cmd/CMakeLists.txt"

#SET(CMD_FILES ${PROJECT_SOURCE_DIR}/cmd/RbGTKGui.cpp)
#ADD_LIBRARY(rb-cmd-lib ${CMD_FILES})

else
echo '
add_executable(rb ${PROJECT_SOURCE_DIR}/revlanguage/main.cpp)

target_link_libraries(rb rb-parser rb-core libs ${Boost_LIBRARIES})
set_target_properties(rb PROPERTIES PREFIX "../")
' >> $HERE/CMakeLists.txt
fi


if [ ! -d "$HERE/libs" ]; then
mkdir "$HERE/libs"
fi
echo 'set(LIBS_FILES' > "$HERE/libs/CMakeLists.txt"
find libs | grep -v "svn" | sed 's|^|${PROJECT_SOURCE_DIR}/|g' >> "$HERE/libs/CMakeLists.txt"
echo ')
add_library(libs ${LIBS_FILES})'  >> "$HERE/libs/CMakeLists.txt"

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

