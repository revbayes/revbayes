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
c = ['ccache', 'x86_64-w64-mingw32-gcc-posix']
cpp = ['ccache', 'x86_64-w64-mingw32-g++-posix']
EOF
)
else
    COMPILER_LINES=$(cat <<EOF
c = 'x86_64-w64-mingw32-gcc-posix'
cpp = 'x86_64-w64-mingw32-g++-posix'
EOF
)
fi

cat > "${CROSSNAME}" <<EOF
[binaries]
${COMPILER_LINES}
ar = 'x86_64-w64-mingw32-ar'
strip = 'x86_64-w64-mingw32-strip'
pkg-config = 'pkg-config'
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

# libwinpthread-git-9.0.0.6090.ad98746a-1
# We can use the libwinpthread-1 from the cross-compiler instead of download it.

PKGS="boost-1.88.0-1
boost-libs-1.88.0-1
dlfcn-1.4.2-1
"

if [ "${gtk}" = "true" ] ; then
    PKGS="$PKGS
    atk-2.56.2-1
    brotli-1.1.0-4
    bzip2-1.0.8-3
    cairo-1.18.4-1
    expat-2.7.1-2
    fontconfig-2.16.2-1
    freetype-2.13.3-1
    fribidi-1.0.16-1
    gdk-pixbuf2-2.42.12-4
    gettext-0.22.4-3
    glib2-2.84.1-2
    graphite2-1.3.14-3
    gtk2-2.24.33-6
    harfbuzz-11.2.0-1
    iconv-1.18-1
    libiconv-1.18-1
    jasper-4.2.5-2
    jbigkit-2.1-5
    lerc-4.0.0-1
    libdatrie-0.2.13-3
    libdeflate-1.23-1
    libffi-3.4.8-1
    libjpeg-turbo-3.1.0-1
    libpng-1.6.48-1
    libthai-0.1.29-3
    libtiff-4.7.0-1
    libwebp-1.5.0-1
    lzo2-2.10-2
    pango-1.56.3-1
    pcre2-10.45-1
    pixman-0.46.0-1
    xz-5.8.1-1
    zlib-1.3.1-1
    zstd-1.5.7-1
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
