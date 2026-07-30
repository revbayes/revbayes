# Building with meson

[Meson](https://mesonbuild.com/) ([git](https://github.com/mesonbuild/meson)) is a build system that fills the same role as CMake.

Meson makes it easy to cross-compile to Windows, Macintosh, and ARM systems from Linux (see [Cross-Compiling](#Cross-compiling))

## Install tools

Before compiling revbayes, you need to install a C++ compiler.  Then you need to install meson and BOOST.

### Mac:
First install XCode and homebrew.

Then using homebrew to install meson and boost:

    brew install meson boost

### Debian/Ubuntu Linux

    apt-get install g++ meson libboost-all-dev


### Redhat Linux

    dnf install gcc-c++ meson boost-devel


### Other systems or not root - use python3

1. You need pip to install meson. Either pip or pip3 is fine, as long as it's using python 3.7 or later.

        pip -V
        pip3 -V

2. Install meson and ninja

        pip3 install --user meson ninja

3. Make meson visible to your shell:

        export PATH=~/.local/bin:$PATH


## Install BOOST

If you do not have superuser access you might not be able to install BOOST using a package manager.
If you are on a computing cluster, there may be a BOOST "module" that you can activate.
If the BOOST library is not installed in a system directory, the `BOOST_ROOT` environment variable will be used to locate it.
To check if this environment variable is set, do:

    echo BOOST_ROOT=${BOOST_ROOT}


### Compile BOOST
If you do not have BOOST installed, you will need to compile it yourself:

    curl -O -L https://archives.boost.io/release/1.88.0/source/boost_1_88_0.tar.gz
    tar -xzvf boost_1_88_0.tar.gz
    cd boost_1_88_0
    ./bootstrap.sh --with-libraries=atomic,chrono,filesystem,system,regex,thread,date_time,program_options,math,serialization --prefix=../installed-boost-1.88.0
    ./b2 link=static cxxflags=-std=c++17 install
    export BOOST_ROOT="$(cd ../installed-boost-1.88.0; pwd)"

The last line above sets `BOOST_ROOT` so that meson can find it.

## Build RevBayes

1. Download RevBayes from our github repository. Clone the repository using git by running the following commands in the terminal:

        git clone https://github.com/revbayes/revbayes.git revbayes
        cd revbayes
        git checkout development           # Probably you want the development branch

2. Configure and compile Revbayes.

    If you are not using the system version of boost, first [tell meson where to find BOOST](#Telling-meson-where-to-find-BOOST).

        projects/meson/generate.sh
        meson setup build --prefix=$HOME/Applications/revbayes
        ninja -C build install

    This creates a `build` directory where the build will take place, and says that the executables will eventually be installed in `$HOME/Applications/revbayes`.

3. For the MPI version, add the `-Dmpi=true` option to the `meson` command.

        meson setup build -Dmpi=true --prefix=$HOME/Applications/revbayes

See the [Configuration options](#Configuration-options) section for more options.


### Troubleshooting
See [If meson cannot find boost](#If-meson-cannot-find-boost) below if this is not working.


## Configuration options

Settings such as `mpi` or `prefix` can be changed using `meson configure`:

    cd build
    meson configure -Dmpi=false
    meson configure -Dprefix=/usr/local

You can also examine current options by leaving off the argument:

    cd build
    meson configure


### Option: MPI

To build with MPI, set the `mpi` option:

    meson setup build-mpi -Dmpi=true -Dprefix=$HOME/Applications/revbayes-mpi
    ninja -C build-mpi install


### Option: RevStudio

To build the `RevStudio` GUI in addition to `rb`, you need to
* install the GTK2 library
* add the `-Dstudio=true` flag when running meson:

To install GTK2 on Debian & Ubuntu:

    apt-get install libgtk+-2.0-dev


To install GTK2 on Mac:

    brew install gtk2


After GTK2 is installed, you can then configure and build RevStudio as follows:

    meson setup build-gtk -Dstudio=true -Dprefix=$HOME/Applications/revbayes-gui
    ninja -C build-gtk install

## Cross-compiling

1. Install the cross-compiler MINGW, if you don't have it already:

    For Debian/Ubuntu Linux:

        apt-get install g++-mingw-w64 wine64-development

    The name of the compiler that we want to use is `x86_64-w64-mingw32-g++-posix`.
    * The prefix `x86_64-w64-` indicates that we want to compile for 64-bit windows, not 32-bit windows.
    * The ending `-posix` indicates that we want to use UNIX (POSIX) threads, not Windows threads.

    There is also a homebrew package `mingw-w64` on Mac, but these instructions have not been tested on Mac.

        brew install mingw-w64 wine-stable

2. Create a fake windows root directory and download windows libraries from MINGW64.

        cd projects/meson/
        ./make_winroot.sh -ccache true

   This command should also generate a "cross-file" that contains info for cross-compiling to windows.

   If you want to compile the GUI version for windows, add `-gtk true` to the command to `make_winroot.sh`.

3. Compile Revbayes

   In the `projects/meson` directory, run:

        ./generate.sh
        meson setup build ../../ --prefix=${HOME}/winrb  --cross-file=win64-cross.txt
        ninja -C build install

   The binary will end up in `~/winrb/bin`.

   If you want to compile the GUI version for windows, add `-Dstudio=true` to the meson command line.

4. Copy needed DLLs to the same directory as `rb.exe`:


        CXX=x86_64-w64-mingw32-g++-posix
        cp -n ~/win_root/mingw64/bin/*.dll                  ~/winrb/bin
        cp -n $($CXX --print-file-name libgcc_s_seh-1.dll)  ~/winrb/bin
        cp -n $($CXX --print-file-name libstdc++-6.dll)     ~/winrb/bin
        cp -n $($CXX --print-file-name libssp-0.dll)        ~/winrb/bin
        cp -n $($CXX --print-file-name libwinpthread-1.dll) ~/winrb/bin

   Note that this same procedure is shown in the `build.yml` workflow.

   If you run `rb.exe` and it cannot find a DLL, it will tell you the first one that it cannot find.


## Troubleshooting:

### If you see `rb: command not found`

The problem is that you tried to run RevBayes but your computer doesn't know where the executable is. The easiest way is to add the directory in which you _installed_ RevBayes to your system path:


    export PATH=<prefix>/bin:$PATH  

For example, if you set `-Dprefix=$HOME/Applications/revbayes`, you would do

    export PATH=$HOME/Applications/revbayes/bin:$PATH

### If meson fails to configure

You can look in `build/meson-logs/meson-log.txt` error log messages to help diagnose the problem.

### If meson cannot find GTK2

If GTK2 is not installed or not usable, this will produce an error message.

Meson uses `pkg-config` to find most dependencies, and to find out what compiler and linker flags they need.

If meson cannot find GTK2, you can:
* run `pkg-config gtk+-2.0 --cflags` to see if pkg-config can find GTK2.
* look in `build-gtk/meson-logs/meson-log.txt` for log messages to help diagnose the problem.

### If meson cannot find boost

First, if `BOOST_ROOT` is set, then check that is correctly describes where the BOOST files are:
  * Run `ls "${BOOST_ROOT}"/include/boost/`.  You should see the boost include files.
  * Run `ls "${BOOST_ROOT}"/lib/`.  You should see the boost library files.

If `ls` shows the required files, then that is a good sign.
However, if the directory does not exist, or if the directory is empty, then you have specified
the environment variables incorrectly.

Second, look in `build/meson-logs/meson-log.txt` for log messages about BOOST.
These will indicate what problem meson is running into.
