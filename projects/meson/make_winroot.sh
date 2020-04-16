#!/bin/bash

# parse command line arguments
while echo $1 | grep ^- > /dev/null; do
    # intercept help while parsing "-key value" pairs
    if [ "$1" = "--help" ] || [ "$1" = "-h" ]
    then
        echo 'Command line options are:
-h                              : print this help and exit.
-ccache          <true|false>   : set to true to add ccache to crossfile

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

# 2. Generate pkg-config wrapper
echo
echo "2. Generating pkg-config wrapper"
cat > "${SYSROOT}/bin/pkg-config" <<EOF
#!/bin/bash
export PKG_CONFIG_SYSROOT_DIR=${SYSROOT}
export PKG_CONFIG_LIBDIR=${SYSROOT}/mingw64/lib/pkgconfig

eval pkg-config "\$@"
EOF
chmod +x "${SYSROOT}/bin/pkg-config"

# 3. Generate cross file
CROSSNAME=win64-cross.txt
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

echo
echo "3. Writing cross file to '${CROSSNAME}'"
cat > "${CROSSNAME}" <<EOF
[binaries]
${COMPILER_LINES}
ar = '/usr/bin/x86_64-w64-mingw32-ar'
strip = '/usr/bin/x86_64-w64-mingw32-strip'
pkgconfig = '${SYSROOT}/bin/pkg-config'
exe_wrapper = 'wine64' # A command used to run generated executables.

# We need these compiler args to find BOOST, which doesn't use pkg-config
[properties]
c_args = ['-I${SYSROOT}/mingw64/include']
c_link_args = ['-L${SYSROOT}/mingw64/lib']

cpp_args = ['-I${SYSROOT}/mingw64/include']
cpp_link_args = ['-L${SYSROOT}/mingw64/lib']
 

[host_machine]
system = 'windows'
cpu_family = 'x86_64'
cpu = 'x86_64'
endian = 'little'
EOF

# 4. Download boost
echo
echo "4. Installing boost to ${SYSROOT}"
echo
cd ${SYSROOT}
PKGS="boost-1.70.0-2"
for PKG in ${PKGS} ; do
    FILE=mingw-w64-x86_64-${PKG}-any.pkg.tar.xz
    rm -f ${FILE}
    wget --no-verbose --show-progress http://repo.msys2.org/mingw/x86_64/${FILE}
    tar -Jxf ${FILE}
    rm ${FILE}
done

echo
echo "Done."
echo
