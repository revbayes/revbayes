# Adapted from https://gist.github.com/peterspackman/8cf73f7f12ba270aa8192d6911972fe8

set(CMAKE_SYSTEM_NAME Windows)
set(TOOLCHAIN_PREFIX x86_64-w64-mingw32)

set(CMAKE_C_COMPILER ${TOOLCHAIN_PREFIX}-gcc)
set(CMAKE_CXX_COMPILER ${TOOLCHAIN_PREFIX}-g++)
set(CMAKE_RC_COMPILER ${TOOLCHAIN_PREFIX}-windres)

# without this you get errors about that RbRegister_VectorType has 'too many sections (33681)'
set(CMAKE_CXX_FLAGS "-Wa,-mbig-obj")
# location of static boost libraries intalled via cygwin
set(Boost_INCLUDE_DIR "/usr/${TOOLCHAIN_PREFIX}/sys-root/mingw/include/")
set(Boost_LIBRARY_DIRS "/usr/${TOOLCHAIN_PREFIX}/sys-root/mingw/lib/")
# required for pkgconfig to successfully find gtk+2.0
set(ENV{PKG_CONFIG_PATH} "/usr/${TOOLCHAIN_PREFIX}/sys-root/mingw/lib/pkgconfig")

# target environment on the build host system
set(CMAKE_FIND_ROOT_PATH "/usr/${TOOLCHAIN_PREFIX}")

# modify default behavior of FIND_XXX() commands
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)