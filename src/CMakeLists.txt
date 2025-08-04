cmake_minimum_required(VERSION 3.5.1)
project(RevBayes)

# Consider:
# -Wno-sign-compare
# -D_GLIBCXX_DEBUG

set(CMAKE_CXX_STANDARD 17)



# Policies relevant to the configuration of Boost libraries
if(POLICY CMP0111)
  # Provides clear error if imported library not found
  cmake_policy(SET CMP0111 NEW)
endif()
if(POLICY CMP0144)
  # Recognize BOOST_ROOT as well as Boost_ROOT
  cmake_policy(SET CMP0144 NEW)
endif()
if(POLICY CMP0167)
  # FindBoost module not available in cmake 3.30
  cmake_policy(SET CMP0167 NEW)
endif()

if(WIN32)
  # Congifure compiler to support big objects
  if(MSVC)
    add_compile_options(/W4 /EHsc /bigobj)
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-mbig-obj")
  endif()
endif()

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

string(TOUPPER "${CMAKE_BUILD_TYPE}" BUILD_TYPE_UPPER)
if (BUILD_TYPE_UPPER STREQUAL "DEBUG")
   set(OPT_FLAGS "${OPT_FLAGS} -g -O0 -Wall")

   # These just add noise:
   set(OPT_FLAGS "${OPT_FLAGS} -Wno-sign-compare -Wno-unknown-pragmas")
   # These could be bugs, but there are so many we can't see the really bad ones. So>
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-reorder")
   # Not necessarily a bug:
   set(OPT_FLAGS "${OPT_FLAGS} -Wno-unused-variable")
# elseif ("${CONTINUOUS_INTEGRATION}" STREQUAL "TRUE")
#    set(OPT_FLAGS "${OPT_FLAGS} -O2 -ffp-contract=off")
#    if (CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|x64|i686|x32|AMD64")
#      set(OPT_FLAGS "${OPT_FLAGS} -mfpmath=sse")
#   endif()
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

## Prevent the error "Could NOT find Threads" when we search for BOOST.
IF(APPLE)
    set(CMAKE_THREAD_LIBS_INIT "-lpthread")
    set(CMAKE_HAVE_THREADS_LIBRARY 1)
    set(CMAKE_USE_WIN32_THREADS_INIT 0)
    set(CMAKE_USE_PTHREADS_INIT 1)
    set(THREADS_PREFER_PTHREAD_FLAG ON)
ENDIF()

##### rpath: where to find shared libraries #####

if ("${CMAKE_SYSTEM_NAME}" MATCHES "Linux")
  # This allows an installed structure of
  #
  # RevBayes/
  #    bin/
  #       rb
  #    lib/
  #       RevBayes/
  #          *.so
  # To find the shared libs from inside bin/ we use $ORIGIN/../lib/RevBayes/
  #
  set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib/RevBayes/")

elseif ("${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")
  # If ever we make an App bundle, the structure might be
  #
  # RevBayes.app/
  #    Contents/
  #       MacOS/
  #          rb
  #       PlugIns/
  #          *.dylib
  #
  # So, at that point, we might add "@loaderpath/../Plugins/"
  #
  set(CMAKE_INSTALL_RPATH "@loaderpath/../lib/RevBayes/;@loaderpath/../PlugIns/")

endif ()

# When user specifies boost paths to use, we need to disable searching from the
# environment and from using cmake's internal boost variables.
if (BOOST_ROOT OR BOOST_INCLUDEDIR OR BOOST_LIBRARYDIR OR DEFINED ENV{BOOST_ROOT} OR DEFINED ENV{BOOST_INCLUDEDIR} OR DEFINED ENV{BOOST_LIBRARYDIR})
  # Don't search for system include dirs for boost
  set(Boost_NO_SYSTEM_PATHS ON)
  # Disable excessive warning messages when using multiple boost versions
  set(Boost_NO_WARN_NEW_VERSION ON)
endif()

# Use the OLD way of looking for Boost
# The NEW way only respects BOOST_ROOT, not these more-specific variables.
if (BOOST_INCLUDEDIR OR BOOST_LIBRARYDIR OR DEFINED ENV{BOOST_INCLUDEDIR} OR DEFINED ENV{BOOST_LIBRARYDIR})
  # Ignore cmake boost variables
  set(Boost_NO_BOOST_CMAKE ON)
endif()

# Report the options that the control how we are looking for boost
MESSAGE("Boost location hints:")
if (BOOST_ROOT)
  MESSAGE("  BOOST_ROOT: ${BOOST_ROOT}")
endif()
if (DEFINED ENV{BOOST_ROOT})
  MESSAGE("  ENV{BOOST_ROOT}: $ENV{BOOST_ROOT}")
endif()
if (NOT DEFINED BOOST_ROOT AND NOT DEFINED ENV{BOOST_ROOT})
  MESSAGE("  BOOST_ROOT not set")
endif()
if (BOOST_INCLUDEDIR)
  MESSAGE("  BOOST_INCLUDEDIR: ${BOOST_INCLUDEDIR}")
endif()
if (DEFINED ENV{BOOST_INCLUDEDIR})
  MESSAGE("  ENV{BOOST_INCLUDEDIR}: $ENV{BOOST_INCLUDEDIR}")
endif()
if (BOOST_LIBRARYDIR)
  MESSAGE("  BOOST_LIBRARYDIR: ${BOOST_LIBRARYDIR}")
endif()
if (DEFINED ENV{BOOST_LIBRARYDIR})
  MESSAGE("  ENV{BOOST_LIBRARYDIR}: $ENV{BOOST_LIBRARYDIR}")
endif()

set(Boost_USE_MULTITHREADED ON)
if ("${STATIC_BOOST}" STREQUAL "ON")
	MESSAGE("Static linking Boost")
	set(Boost_USE_STATIC_LIBS ON)
endif()

message("Searching for BOOST:")
find_package(Boost
1.74.0
COMPONENTS regex
program_options
thread
system
date_time
serialization REQUIRED)
message("Boost configuration results:")
message("  Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
message("  Boost_LIBRARY_DIRS: ${Boost_LIBRARY_DIRS}")
message("  Boost_LIBRARIES: ${Boost_LIBRARIES}")
link_directories(${Boost_LIBRARY_DIRS})

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

# Should we build help2yml?
if ("${HELP}" STREQUAL "ON")
  add_subdirectory(help2yml)

  message("Building ${RB_EXEC_NAME}-help2yml")
  add_executable(${RB_EXEC_NAME}-help2yml ${PROJECT_SOURCE_DIR}/help2yml/main.cpp)

  target_link_libraries(${RB_EXEC_NAME}-help2yml rb-help rb-parser rb-core rb-libs rb-parser ${Boost_LIBRARIES} ${CMAKE_DL_LIBS})
  set_target_properties(${RB_EXEC_NAME}-help2yml PROPERTIES PREFIX "../")
  if ("${MPI}" STREQUAL "ON")
    target_link_libraries(${RB_EXEC_NAME}-help2yml ${MPI_LIBRARIES})
  endif()
endif()

if ("${CMD_GTK}" STREQUAL "ON")
  set(RB_EXEC_NAME "RevStudio")
endif()
message("Building ${RB_EXEC_NAME}")

if ("${CMD_GTK}" STREQUAL "ON")
  target_link_libraries(rb-jupyter rb-parser rb-core rb-libs ${Boost_LIBRARIES} ${CMAKE_DL_LIBS})

  # Use the package PkgConfig to detect GTK+ headers/library files
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(GTK REQUIRED gtk+-2.0)
  #pkg_check_modules(GTK REQUIRED gtk+-3.0)

  # Setup CMake to use GTK+, tell the compiler where to look for headers
  # and to the linker where to look for libraries
  include_directories(${GTK_INCLUDE_DIRS})
  link_directories(${GTK_LIBRARY_DIRS})
  add_definitions(${GTK_CFLAGS_OTHER})      # Add other flags to the compiler
  TARGET_LINK_LIBRARIES(RevStudio rb-cmd-lib rb-parser rb-core rb-libs ${Boost_LIBRARIES} ${GTK_LIBRARIES} ${CMAKE_DL_LIBS})

  # Look into the cmd directory
  add_subdirectory(cmd)
  set(RB_CMD_LIB rb-cmd-lib)

  # The GTK version uses a different main.cpp
  add_executable(${RB_EXEC_NAME}   ${PROJECT_SOURCE_DIR}/cmd/main.cpp)

else()
  add_executable(${RB_EXEC_NAME}   ${PROJECT_SOURCE_DIR}/revlanguage/main.cpp)
endif()

# Link the target to other libraries.  Some variables can be empty.
target_link_libraries(${RB_EXEC_NAME}
  ${RB_CMD_LIB}
  rb-parser
  rb-core
  rb-libs
  ${Boost_LIBRARIES}
  ${CMAKE_DL_LIBS}
  ${GTK_LIBRARIES}
  ${MPI_LIBRARIES}
)

install(TARGETS ${RB_EXEC_NAME} DESTINATION bin)
