# Building with meson

[Meson](https://mesonbuild.com/) ([git](https://github.com/mesonbuild/meson)) is a new build system that aims to replace autotools and fills the same role as CMake.

Meson makes it easy to cross-compile to Windows, Macintosh, and ARM systems from Linux (see [Cross-Compiling](#Cross-compiling))

## Install tools

Before compiling revbayes, you need to install a C++ compiler.  Then you need to install meson and BOOST.

### Mac:
- First install XCode and homebrew
``` sh
brew install meson boost
```

### Debian/Ubuntu Linux
``` sh
apt-get install g++ meson libboost-all-dev
```

### Redhat Linux
``` sh
dnf install gcc-c++ meson boost-devel
```

### Other systems or not root - use python3

1. You need pip to install meson. Either pip or pip3 is fine, as long as its using python 3.7 or later.

   ``` sh
   pip -V
   pip3 -V
   ```

2. Install meson and ninja

   ``` sh
   pip3 install --user meson ninja
   ```

3. Make meson visible to your shell:

   ``` sh
   export PATH=~/.local/bin:$PATH
   ```

## Build RevBayes

1. Download RevBayes from our github repository. Clone the repository using git by running the following commands in the terminal:

    ``` sh
    git clone https://github.com/revbayes/revbayes.git revbayes
    cd revbayes
    git checkout development           # Probably you want the development branch
    git submodule init
    git submodule update               # Download tests required for the build
    ```

2. Configure and compile Revbayes.

    If you are not using the system version of boost, first [tell meson where to find BOOST](#Telling-meson-where-to-find-BOOST).
   
    ``` sh
    projects/meson/generate.sh
    meson build --prefix=$HOME/Applications/revbayes
    ninja -C build install
    ```

    This creates a `build` directory where the build will take place, and says that the executables will eventually be installed in `$HOME/Applications/revbayes`.

3. For the MPI version, add the `-Dmpi=true` option to the `meson` command.

    ``` sh
    meson build -Dmpi=true --prefix=$HOME/Applications/revbayes
    ```

    See the [Configuration options](#Configuration-options) section for more options.

## Telling meson where to find BOOST

If you need to use a version of BOOST that is not the system version, you can specify the location of BOOST 
by setting environment variables before you run meson.

### Checking BOOST environment variables
If you are on a cluster, these variables might be already set.
Or they might get set if you load a boost "module".
To check if they are set and what their values are, run
``` sh
echo BOOST_ROOT=${BOOST_ROOT}
echo BOOST_INCLUDEDIR=${BOOST_INCLUDEDIR}
echo BOOST_LIBRARYDIR=${BOOST_LIBRARYDIR}
```
You can also use the `unset` command to unset them.  For example: `unset BOOST_ROOT`.

### Setting BOOST environment variables -- simple
Many simple installations have include files in `/path/to/boost/include`, and library files in
`/path/to/boost/lib`.  In this case you can specify:
``` sh
export BOOST_ROOT=/path/to/boost/
 ```

### Setting BOOST environment variables -- general
In some cases you need to specify the directories that contain include files and library files directly.
For example, suppose that you have compiled your own boost in `/path/to/boost_1_74_0`.
In this case, you can specify
``` sh
export BOOST_INCLUDEDIR=/path/to/boost_1_74_0/
export BOOST_LIBRARYDIR=/path/to/boost_1_74_0/stage/lib
```
These variables override BOOST_ROOT.

### Troubleshooting
See [If meson cannot find boost](#If-meson-cannot-find-boost) below if this is not working.

## Configuration options

Settings such as `mpi` or `prefix` can be changed using `meson configure`:
``` sh
cd build
meson configure -Dmpi=false
meson configure -Dprefix=/usr/local
```
You can also examine current options by leaving off the argument:
``` sh
cd build
meson configure
```

### Option: MPI

To build with MPI, set the `mpi` option:
```
meson build-mpi -Dmpi=true -Dprefix=$HOME/Applications/revbayes-mpi
ninja -C build-mpi install
```

### Option: Jupyter Kernel

To build the jupyter kernel, set the `jupyter` option:
```
meson build-jupyter -Djupyter=true -Dprefix=$HOME/Applications/revbayes-jupyter
ninja -C build-jupyter install
```

### Option: RevStudio

To build the `RevStudio` GUI in addition to `rb`, you need to
* install the GTK2 library
* add the `-Dstudio=true` flag when running meson:

To install GTK2 on Debian & Ubuntu:
```
apt-get install libgtk+-2.0-dev
```

To install GTK2 on Mac:
```
brew install gtk2
```

After GTK2 is installed, you can then configure and build RevStudio as follows:
```
meson build-gtk -Dstudio=true -Dprefix=$HOME/Applications/revbayes-gui
ninja -C build-gtk install
```

## Cross-compiling

1. Install the cross-compiler MINGW, if you don't have it already:

    For Debian/Ubuntu Linux:

        apt-get install g++-mingw-w64 wine64-development

    The name of the compiler that we want to use is `x86_64-w64-mingw32-g++-posix`.
    * The prefix `x86_64-w64-` indicates that we want to compile for 64-bit windows, not 32-bit windows.
    * The ending `-posix` indicates that we want to use UNIX (POSIX) threads, not Windows threads.

    There is also a homebrew package `mingw-w64` on Mac, but these instructions have not been tested on Mac.

        brew install mingw-w64 wine-stable

    You might need to update the paths in the cross-file below.

2. Create a fake windows root directory and download windows libraries from MINGW64.

   ```
   cd projects/meson/
   ./make_winroot.sh -ccache true
   ```

   This command should also generate a "cross-file" that contains info for cross-compiling to windows.

   If you want to compile the GUI version for windows, add `-gtk true` to the command to `make_winroot.sh`.


3. Compile Revbayes

   In the `projects/meson` directory, run:

   ``` sh
   ./generate.sh
   meson build ../../ --prefix=${HOME}/winrb  --cross-file=win64-cross.txt
   ninja -C build install
   ```

   The binary will end up in `~/winrb/bin`.

   If you want to compile the GUI version for windows, add `-Dstudio=true` to the meson command line.

4. Copy needed DLLs to the same directory as `rb.exe`:

   ```
   CXX=x86_64-w64-mingw32-g++-posix
   cp -n ~/win_root/mingw64/bin/*.dll                  ~/winrb/bin
   cp -n $($CXX --print-file-name libgcc_s_seh-1.dll)  ~/winrb/bin
   cp -n $($CXX --print-file-name libstdc++-6.dll)     ~/winrb/bin
   cp -n $($CXX --print-file-name libssp-0.dll)        ~/winrb/bin
   cp -n $($CXX --print-file-name libwinpthread-1.dll) ~/winrb/bin
   ```

   Note that this same procedure is shown in the `build.yml` workflow.

   If you run `rb.exe` and it cannot find a DLL, it will tell you the first one that it cannot find.



## Troubleshooting:

### If you see `rb: command not found`

The problem is that you tried to run RevBayes but your computer doesn't know where the executable is. The easiest way is to add the directory in which you _installed_ RevBayes to your system path:

```
export PATH=<prefix>/bin:$PATH  
```

For example, if you set `-Dprefix=$HOME/Applications/revbayes`, you would do

```
export PATH=$HOME/Applications/revbayes/bin:$PATH
```

### If meson fails to configure

You can look in `build/meson-logs/meson-log.txt` error log messages to help diagnose the problem.

### If meson cannot find GTK2

If GTK2 is not installed or not usable, this will produce an error message.

Meson uses `pkg-config` to find most dependencies, and to find out what compiler and linker flags they need.

If meson cannot find GTK2, you can:
* run `pkg-config gtk+-2.0 --cflags` to see if pkg-config can find GTK2.
* look in `build-gtk/meson-logs/meson-log.txt` for log messages to help diagnose the problem.

### If meson cannot find boost

First, check that you correctly told meson where to find the boost files.
* if you specified BOOST_ROOT, then
    * run `ls ${BOOST_ROOT}/include/boost/`.  You should see the boost include files.
    * run `ls ${BOOST_ROOT}/lib/`.  You should see boost library files.
* if you specified BOOST_INCLUDEDIR and BOOST_LIBRARYDIR, then
    * run `ls ${BOOST_INCLUDEDIR}/boost/`.  You should see the boost include files.
    * run `ls ${BOOST_LIBRARYDIR}`.  You should see boost library files.

If `ls` shows the required files, then that is a good sign.
However, if the directory does not exist, or if the directory is empty, then you have specified
the environment variables incorrectly.

Second, look in `build/meson-logs/meson-log.txt` for log messages to help diagnose the problem.
They should be in a section that begins with `Trying to find boost with:`.

A successful run with BOOST_INCLUDEDIR and BOOST_LIBRARYDIR specified looks something like:
```
Trying to find boost with:
  - boost_includedir = /path/to/boost_1_74_0
  - boost_librarydir = /path/to/boost_1_74_0/stage/lib
  - potential library dirs: ['/path/to/boost_1_74_0/stage/lib']
  - potential include dirs: ['/path/to/boost_1_74_0']
  - found boost library dir: /path/to/boost_1_74_0/stage/lib
  - found boost 1.74.0 include dir: /path/to/boost_1_74_0
  - filtered library list:
    - <LIB: -M ------ ??? ? 1_74 boost_atomic                     /path/to/boost_1_74_0/stage/lib/libboost_atomic.so.1.74.0>
  ...
```
