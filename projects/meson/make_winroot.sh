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
-gtk             <true|false>   : include GTK2 packages.

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

if [ "$ccache" != "false" ]; then
    COMPILER_LINES=$(cat <<EOF
c = ['ccache', '/usr/bin/x86_64-w64-mingw32-gcc-posix']
cpp = ['ccache', '/usr/bin/x86_64-w64-mingw32-g++-posix']
EOF
)
else
    COMPILER_LINES=$(cat <<EOF
c = '/usr/bin/x86_64-w64-mingw32-gcc-posix'
cpp = '/usr/bin/x86_64-w64-mingw32-g++-posix'
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

# Note that the use of gcc-posix and g++-posix means that we need
# *-posix/libgcc_s_seh-1.dll and *-posix2/libstdc++-6.dll instead
# of the *-win32/ versions.

PKGS="boost-1.75.0-2
libwinpthread-git-9.0.0.6090.ad98746a-1
openlibm-0.7.5-1
"

if [ "${gtk}" = "true" ] ; then
    PKGS="$PKGS
    atk-2.36.0-2
    brotli-1.0.9-4
    bzip2-1.0.8-2
    cairo-1.16.0-3
    expat-2.4.6-1
    fontconfig-2.13.96-1
    freetype-2.11.1-2
    fribidi-1.0.11-1
    gdk-pixbuf2-2.42.6-2
    gettext-0.21-3
    glib2-2.70.4-1
    graphite2-1.3.14-2
    gtk2-2.24.33-4
    harfbuzz-3.4.0-1
    iconv-1.16-2
    jasper-2.0.33-1
    libdatrie-0.2.13-1
    libffi-3.3-4
    libiconv-1.16-2
    libjpeg-turbo-2.1.3-1
    libpng-1.6.37-6
    libthai-0.1.29-1
    libtiff-4.3.0-7
    pango-1.50.4-1
    pcre-8.45-1
    pixman-0.40.0-2
    xz-5.2.5-2
    zlib-1.2.11-9
"
fi

for PKG in ${PKGS} ; do
    FILE=mingw-w64-x86_64-${PKG}-any.pkg.tar.zst
    if [ -e "${FILE}" ] ; then
        echo "   ${PKG} already downloaded and installed."
    else
        URL="http://repo.msys2.org/mingw/x86_64/${FILE}"
        if ! wget --no-verbose --show-progress "${URL}" ; then
            echo "Failed to download ${URL}"
            exit
        fi
        if tar -I zstd -xf ${FILE} ; then
            echo "${PKG} installed"
        else
            rm ${FILE}
            echo "Failed to install ${PKG}"
            echo
            echo "Perhaps the decompression program zstd is not installed?"
            echo "Try 'sudo apt install zstd'."
            exit
        fi
    fi
done
fi

echo
echo "Done."
echo
