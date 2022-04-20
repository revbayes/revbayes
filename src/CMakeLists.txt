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

if ("${JUPYTER}" STREQUAL "ON")
   add_definitions(-DRB_XCODE)
endif()

# When user specifies boost paths to use, we need to disable searching from the
# environment and from using cmake's internal boost variables.
if (BOOST_ROOT OR BOOST_INCLUDEDIR OR BOOST_LIBRARYDIR) 
  # Ignore environment boost variables
  set(Boost_NO_SYSTEM_PATHS ON)
  # Ignore cmake boost variables
  set(Boost_NO_BOOST_CMAKE ON)
  # Disable excessive warning messages when using multiple boost versions
  set(Boost_NO_WARN_NEW_VERSION ON)
endif()

# For standard installations, only `BOOST_ROOT` needs to be set by the user.
if (BOOST_ROOT)
  MESSAGE("  BOOST_ROOT: ${BOOST_ROOT}")
  # We still set the internal boost variables for nicer messages
  set(Boost_INCLUDE_DIR ${BOOST_ROOT}/include)
  set(Boost_LIBRARY_DIR ${BOOST_ROOT}/lib)
endif()

# For non-standard installations where boost's `include` and `lib` directories
# are not located in the same directory, we can specify `BOOST_INCLUDEDIR`
# and # `BOOST_LIBRARYDIR`.
if (BOOST_INCLUDEDIR AND BOOST_LIBRARYDIR)
  # Set `include` directory
  MESSAGE("  BOOST_INCLUDEDIR: ${BOOST_INCLUDEDIR}")
  set(Boost_INCLUDE_DIR ${BOOST_ROOT})
  # Set `lib` directory
  MESSAGE("  BOOST_LIBRARYDIR: ${BOOST_LIBRARYDIR}")
  set(Boost_LIBRARY_DIR ${BOOST_LIBRARYDIR})
# We need both `BOOST_INCLUDEDIR` and `BOOST_LIBRARYDIR` to be specified. 
# These variables should never occur on their own. 
elseif  (BOOST_INCLUDEDIR AND NOT BOOST_LIBRARYDIR)
  MESSAGE("Warning: BOOST_LIBRARYDIR should probably be set.")
elseif  (NOT BOOST_INCLUDEDIR AND BOOST_LIBRARYDIR)
  MESSAGE("Warning: BOOST_INCLUDEDIR should probably be set.")
endif()

# Enable boost multithreading
set(Boost_USE_MULTITHREADED ON)

# Use static linking if user has asked for it.
if ("${STATIC_BOOST}" STREQUAL "ON")
  MESSAGE("Static linking Boost")
  set(Boost_USE_STATIC_LIBS ON)
endif()

# Find required boost packages
MESSAGE("Searching for BOOST:")
find_package(Boost
1.71.0
COMPONENTS regex
program_options
thread
system
filesystem
date_time
serialization REQUIRED)

# Add boost to the directories to include and link.
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

# Print results of the boost configuration
MESSAGE("Boost configuration results:")
MESSAGE("  Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIR}")
MESSAGE("  Boost_LIBRARY_DIRS: ${Boost_LIBRARY_DIRS}")
MESSAGE("  Boost_LIBRARIES: ${Boost_LIBRARIES}")


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

