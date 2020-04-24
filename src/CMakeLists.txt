cmake_minimum_required(VERSION 2.6)
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

# Add extra CMake libraries into ./CMake
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR} ${CMAKE_MODULE_PATH})

set(OPTIMIZATION_FLAGS "-O3 -DNDEBUG")
if (CMAKE_SYSTEM_PROCESSOR MATCHES "^arm*|aarch64")
   add_definitions(-DRB_ARM)
else()
   set(OPTIMIZATION_FLAGS "${OPTIMIZATION_FLAGS} -msse -msse2 -msse3")
endif()

if (CMAKE_BUILD_TYPE MATCHES DEBUG)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -Wall")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O0 -Wall")
else()
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPTIMIZATION_FLAGS}")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OPTIMIZATION_FLAGS}")
endif()

if ("${RB_EXEC_NAME}" STREQUAL "")
   SET(RB_EXEC_NAME "rb")
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

if (WIN32)
   add_definitions(-DRB_WIN)
endif()

if ("${JUPYTER}" STREQUAL "ON")
   add_definitions(-DRB_XCODE)
endif()

if ("${CONTINUOUS_INTEGRATION}" STREQUAL "TRUE")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g0 -O2")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g0 -O2")
endif()


MESSAGE("My Boost information:")
MESSAGE("  Boost_INCLUDE_DIRS: ${LOCAL_BOOST_ROOT}")
MESSAGE("  Boost_LIBRARIES: ${LOCAL_BOOST_LIBRARY}")

set(Boost_USE_MULTITHREADED ON)

if ( NOT "${LOCAL_BOOST_ROOT}" STREQUAL "" AND NOT "${LOCAL_BOOST_LIBRARY}" STREQUAL "" )
   SET(BOOST_ROOT ${LOCAL_BOOST_ROOT})
   SET(BOOST_LIBRARY ${LOCAL_BOOST_LIBRARY})
   SET(Boost_REALPATH ON)
#   SET(Boost_NO_SYSTEM_PATHS ON)
#   SET(Boost_USE_STATIC_RUNTIME ON)
#   SET(Boost_USE_STATIC_LIBS ON)
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

  target_link_libraries(${RB_EXEC_NAME} rb-parser rb-core rb-libs ${Boost_LIBRARIES})

  set_target_properties(${RB_EXEC_NAME} PROPERTIES PREFIX "../")

  if ("${MPI}" STREQUAL "ON")
    target_link_libraries(${RB_EXEC_NAME} ${MPI_LIBRARIES})
  endif()

endif()

