#!/bin/bash

# parse command line arguments
while echo $1 | grep ^- > /dev/null; do
    # intercept help while parsing "-key value" pairs
    if [ "$1" = "--help" ] || [ "$1" = "-h" ]
    then
        echo 'Command line options are:
-h                              : print this help and exit.
-ccache          <true|false>   : set to true to add ccache to crossfile.
-download        <true|false>   : download MINGW packages.

Example:
  ./make_winroot.sh -help true'
        exit
    fi

    # parse pairs
    eval $( echo $1 | sed 's/-//g' | tr -d '\012')=$2
    shift
    shift
done

SYSROOT=$HOME/win_root

# 1. Make sysroot
echo
echo "1. Writing sysroot dir ${SYSROOT}"
mkdir -p "${SYSROOT}"
mkdir -p "${SYSROOT}/bin"

# 2. Generate cross file
CROSSNAME=win64-cross.txt
echo
echo "2. Writing cross file to '${CROSSNAME}'"

if [ "$ccache" = "true" ]; then
    COMPILER_LINES=$(cat <<EOF
c = ['ccache', '/usr/bin/x86_64-w64-mingw32-gcc']
cpp = ['ccache', '/usr/bin/x86_64-w64-mingw32-g++']
EOF
)
else
    COMPILER_LINES=$(cat <<EOF
c = '/usr/bin/x86_64-w64-mingw32-gcc'
cpp = '/usr/bin/x86_64-w64-mingw32-g++'
EOF
)
fi

cat > "${CROSSNAME}" <<EOF
[binaries]
${COMPILER_LINES}
ar = '/usr/bin/x86_64-w64-mingw32-ar'
strip = '/usr/bin/x86_64-w64-mingw32-strip'
pkgconfig = '/usr/bin/pkg-config'
exe_wrapper = 'wine64' # A command used to run generated executables.

# why do we still need these? shouldn't they get added automatically if we find boost?
[built-in options]
c_args = ['-I${SYSROOT}/mingw64/include']
c_link_args = ['-L${SYSROOT}/mingw64/lib']

cpp_args = ['-I${SYSROOT}/mingw64/include']
cpp_link_args = ['-L${SYSROOT}/mingw64/lib']

[properties]
sys_root = '${SYSROOT}'
pkg_config_libdir = '${SYSROOT}/mingw64/lib/pkgconfig'
boost_root = '${SYSROOT}/mingw64'

[host_machine]
system = 'windows'
cpu_family = 'x86_64'
cpu = 'x86_64'
endian = 'little'
EOF

# 3. Download packages
if [ "$download" = "false" ] ; then
    echo
    echo "(3.) Skipping installation of packages to ${SYSROOT}."
    echo
else

    echo
    echo "3. Installing packages to ${SYSROOT}"
    echo
cd ${SYSROOT}
PKGS="boost-1.75.0-2 libwinpthread-git-9.0.0.6090.ad98746a-1 openlibm-0.7.5-1"
for PKG in ${PKGS} ; do
    FILE=mingw-w64-x86_64-${PKG}-any.pkg.tar.zst
    rm -f ${FILE}
    wget --no-verbose --show-progress http://repo.msys2.org/mingw/x86_64/${FILE}
    if tar -I zstd -xf ${FILE} ; then
        rm ${FILE}
    else
        echo
        echo "Perhaps the decompression program zstd is not installed?"
        echo "Try 'sudo apt install zstd'."
        exit
    fi
done
fi

echo
echo "Done."
echo
