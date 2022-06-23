cmake_minimum_required(VERSION 3.5.1)
project(RevBayes)

# Consider:
# -Wno-sign-compare
# -D_GLIBCXX_DEBUG

# This is the RIGHT way, but requires cmake version >=3:
#   set(CMAKE_CXX_STANDARD 11)
#
# https://gitlab.kitware.com/cmake/community/-/wikis/CMake-Versions-on-Linux-Distros
#
# So, we add the flag directly instead.
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
endif()

# Add extra CMake libraries into ./CMake
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR} ${CMAKE_MODULE_PATH})

if (CMAKE_SYSTEM_PROCESSOR MATCHES "^arm*|aarch64")
   add_definitions(-DRB_ARM)

   set(OPT_FLAGS "")
else()
   #  We need to enable SSE3 instructions even in debug builds.
   #  Otherwise we get errors failing to inline __mm_hadd_pd( ) in PhyloCTMCSiteHomogeneousNucleotide.h
   set(OPT_FLAGS "-msse -msse2 -msse3")
endif()

if(CMAKE_BUILD_TYPE MATCHES DEBUG)
   set(OPT_FLAGS "${OPT_FLAGS} -g -O0 -Wall")

   # These just add noise:
   set(OPT_FLAGS "${OPT_FLAGS} -Wno-sign-compare -Wno-unknown-pragmas")
   # These could be bugs, but there are so many we can't see the really bad ones. So>
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-reorder")
   # Not necessarily a bug:
   set(OPT_FLAGS "${OPT_FLAGS} -Wno-unused-variable")
elseif ("${CONTINUOUS_INTEGRATION}" STREQUAL "TRUE")
   set(OPT_FLAGS "${OPT_FLAGS} -O2 -mfpmath=sse -ffp-contract=off")
   find_library(OPENLIBM openlibm)
   message("openlibm: ${OPENLIBM}")
else()
   set(OPT_FLAGS "${OPT_FLAGS} -O3 -DNDEBUG")
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OPT_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPT_FLAGS}")

if ("${RB_EXEC_NAME}" STREQUAL "")
   set(RB_EXEC_NAME "rb")
endif()

if ("${MPI}" STREQUAL "ON")
   add_definitions(-DRB_MPI)
   #add_definitions(-DDEBUG_MPI_MCA)
   # Require MPI for this project:
   find_package(MPI REQUIRED)
   include_directories(${MPI_INCLUDE_PATH})
   set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
   set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS}")
endif()


if ("${RB_BEAGLE}" STREQUAL "ON")
    if (RB_BEAGLE_ROOT)
        set(BEAGLE_ROOT ${RB_BEAGLE_ROOT})
	set(BEAGLE_INCLUDE_DIR, "${BEAGLE_ROOT}/include/libhmsbeagle-1")	
	set(BEAGLE_LIBRARY_DIRS, "${BEAGLE_ROOT}/lib ${BEAGLE_ROOT}/lib/libhmsbeagle-1")	

	set(BEAGLE_COMPILE_FLAGS    "-I${BEAGLE_INCLUDE_DIR}")
	set(BEAGLE_LINK_FLAGS       "-L${BEAGLE_ROOT}/lib/libhmsbeagle-1 -L${BEAGLE_ROOT}/lib -lhmsbeagle")
	set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} ${BEAGLE_COMPILE_FLAGS}")
	set(CMAKE_CXX_LINK_FLAGS    "${CMAKE_CXX_LINK_FLAGS} ${BEAGLE_LINK_FLAGS}")

	#set(BEAGLE_INCLUDE_PATH  "${BEAGLE_ROOT}/include/libhmsbeagle-1")
	#include_directories(${BEAGLE_INCLUDE_PATH})
	#set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} ${BEAGLE_COMPILE_FLAGS}")
	#set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${BEAGLE_LINK_FLAGS}")

	#MESSAGE("My BEAGLE information:")
        #MESSAGE("  BEAGLE_INCLUDE : ${BEAGLE_INCLUDE_PATH}")
        #MESSAGE("  BEAGLE_LIBS    : ${BEAGLE_LINK_FLAGS}")
    else()
	#-- TODO - Check this on other computers...
	MESSAGE("Searching for BEAGLE packages using `pkg-config`:")
        find_package(PkgConfig REQUIRED)
	pkg_check_modules(BEAGLE REQUIRED hmsbeagle-1 hmsbeagle)
    endif()

    add_definitions(-DRB_BEAGLE)
    
    INCLUDE_DIRECTORIES(${BEAGLE_INCLUDE_DIR})
    LINK_DIRECTORIES(${BEAGLE_LIBRARY_DIRS})
	
    MESSAGE("BEAGLE configuration results:")
    MESSAGE("  BEAGLE_INCLUDE_DIR: ${BEAGLE_INCLUDE_DIR}")
    MESSAGE("  BEAGLE_LIBRARY_DIR: ${BEAGLE_LIBRARY_DIRS}")
endif()


if ("${RB_USE_EIGEN3}" STREQUAL "ON")
   find_package(Eigen3 NAMES Eigen3 Eigen eigen3 eigen REQUIRED)
   include_directories(${EIGEN3_INCLUDE_DIR})

   add_definitions(-DRB_USE_EIGEN3)

   set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} ${EIGEN_COMPILE_FLAGS}")
   MESSAGE("My EIGEN3 information:")
   MESSAGE("  EIGEN3_INCLUDE_DIR : ${EIGEN3_INCLUDE_DIR}")
   MESSAGE("  EIGEN3_DIR         : ${EIGEN3}")
endif()


if ("${JUPYTER}" STREQUAL "ON")
   add_definitions(-DRB_XCODE)
endif()

MESSAGE("Boost location hints:")
if (BOOST_ROOT OR BOOST_INCLUDEDIR OR BOOST_LIBRARYDIR) 
  set(Boost_NO_SYSTEM_PATHS ON)
  set(Boost_NO_BOOST_CMAKE ON)
  set(Boost_NO_WARN_NEW_VERSION ON)
endif()
if (BOOST_ROOT)
  MESSAGE("  BOOST_ROOT: ${BOOST_ROOT}")
  set(Boost_INCLUDE_DIR ${BOOST_ROOT}/include)
  set(Boost_LIBRARY_DIR ${BOOST_ROOT}/lib)
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
endif()
if (BOOST_INCLUDEDIR)
  MESSAGE("  BOOST_INCLUDEDIR: ${BOOST_INCLUDEDIR}")
  set(Boost_INCLUDE_DIR ${BOOST_ROOT})
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
endif()
if (BOOST_LIBRARYDIR)
  MESSAGE("  BOOST_LIBRARYDIR: ${BOOST_LIBRARYDIR}")
  set(Boost_LIBRARY_DIR ${BOOST_LIBRARYDIR})
endif()

set(Boost_USE_MULTITHREADED ON)
if ("${STATIC_BOOST}" STREQUAL "ON")
	MESSAGE("Static linking Boost")
	set(Boost_USE_STATIC_LIBS ON)
endif()

MESSAGE("Searching for BOOST packages:")
find_package(Boost
  1.74.0
  COMPONENTS regex
  program_options
  thread
  system
  filesystem
  date_time
  serialization REQUIRED)

MESSAGE("Boost configuration results:")
MESSAGE("  Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIR}")
MESSAGE("  Boost_LIBRARY_DIRS: ${Boost_LIBRARY_DIRS}")
MESSAGE("  Boost_LIBRARIES: ${Boost_LIBRARIES}")
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

# This will look for "generated_include_dirs.cmake" in the module path.
include("generated_include_dirs")

# Split into much smaller libraries
add_subdirectory(libs)
add_subdirectory(core)
add_subdirectory(revlanguage)


message("Compile information:")
message("  CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
message("  CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
message("  CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")

############# executables #################
# basic rev-bayes binary

if ("${HELP}" STREQUAL "ON")
  add_subdirectory(help2yml)

  message("Building ${RB_EXEC_NAME}-help2yml")
  add_executable(${RB_EXEC_NAME}-help2yml ${PROJECT_SOURCE_DIR}/help2yml/main.cpp)

  target_link_libraries(${RB_EXEC_NAME}-help2yml rb-help rb-parser rb-core rb-libs rb-parser ${Boost_LIBRARIES})
  set_target_properties(${RB_EXEC_NAME}-help2yml PROPERTIES PREFIX "../")
  if ("${MPI}" STREQUAL "ON")
    target_link_libraries(${RB_EXEC_NAME}-help2yml ${MPI_LIBRARIES})
  endif()
endif()

if ("${JUPYTER}" STREQUAL "ON")
  message("Building rb-jupyter")
  add_executable(rb-jupyter ${PROJECT_SOURCE_DIR}/revlanguage/main.cpp)

  target_link_libraries(rb-jupyter rb-parser rb-core rb-libs ${Boost_LIBRARIES})
  set_target_properties(rb-jupyter PROPERTIES PREFIX "../")
elseif ("${CMD_GTK}" STREQUAL "ON")
  message("Building RevStudio")
  # Use the package PkgConfig to detect GTK+ headers/library files
  FIND_PACKAGE(PkgConfig REQUIRED)
  PKG_CHECK_MODULES(GTK REQUIRED gtk+-2.0)
  #PKG_CHECK_MODULES(GTK REQUIRED gtk+-3.0)

  # Setup CMake to use GTK+, tell the compiler where to look for headers
  # and to the linker where to look for libraries
  INCLUDE_DIRECTORIES(${GTK_INCLUDE_DIRS})
  LINK_DIRECTORIES(${GTK_LIBRARY_DIRS})

  # Add other flags to the compiler
  ADD_DEFINITIONS(${GTK_CFLAGS_OTHER})

  # Add an executable compiled from hello.c
  ADD_EXECUTABLE(RevStudio ${PROJECT_SOURCE_DIR}/cmd/main.cpp)

  # Link the target to the GTK+ libraries
  TARGET_LINK_LIBRARIES(RevStudio rb-cmd-lib rb-parser rb-core rb-libs ${Boost_LIBRARIES} ${GTK_LIBRARIES})

  SET_TARGET_PROPERTIES(RevStudio PROPERTIES PREFIX "../")

  add_subdirectory(cmd)

else()
  message("Building ${RB_EXEC_NAME}")
  add_executable(${RB_EXEC_NAME} ${PROJECT_SOURCE_DIR}/revlanguage/main.cpp)

  target_link_libraries(${RB_EXEC_NAME} rb-parser rb-core rb-libs ${Boost_LIBRARIES} ${OPENLIBM})

  set_target_properties(${RB_EXEC_NAME} PROPERTIES PREFIX "../")

  if ("${MPI}" STREQUAL "ON")
    target_link_libraries(${RB_EXEC_NAME} ${MPI_LIBRARIES})
  endif()

endif()

