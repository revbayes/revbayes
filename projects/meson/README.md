# Building with meson

[Meson](https://mesonbuild.com/) ([git](https://github.com/mesonbuild/meson)) is a new build system that aims to replace autotools and fills the same role as CMake.

Meson makes it easy to cross-compile to Windows, Macintosh, and ARM systems from Linux (see [Cross-Compile.md](Cross-Compile.md)).

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

### Other

You can install meson using python:
``` sh
python3 -m venv meson_env
source meson_env/bin/activate
pip3 install meson ninja
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

2. Configure and compile Revbayes:
   
    ``` sh
    projects/meson/generate.sh
    meson build --prefix=$HOME/Applications/revbayes
    ninja -C build install
    ```

    This creates a `build` directory where the build will take place, and says that the executables will eventually be installed in `$HOME/Applications/revbayes`.

3. For the MPI version, add `-Dmpi=true` to the `meson` command.

    ``` sh
    meson build -Dmpi=true --prefix=$HOME/Applications/revbayes
    ```

## Configuration

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

## Cross-compiling (Linux -> Windows)

1. Create a fake windows root directory and download windows libraries from MINGW64.

   ```
   cd projects/meson/
   ./make_winroot.sh -ccache true
   ```

   This command should also generate a "cross-file" that contains info for cross-compiling to windows.

   If you want to compile the GUI version for windows, add `-gtk true` to the command to `make_winroot.sh`.

2. Compile Revbayes

   ```
   cd projects/meson
   ./generate.sh
   meson build ../../ --prefix=${HOME}/winrb  --cross-file=win64-cross.txt
   ninja -C build install
   ```

   The binary will end up in `~/winrb/bin`.

   If you want to compile the GUI version for windows, add `-Dstudio=true` to the meson command line.

3. Copy DLLs that the binary needs to the same directory

   ```
   cp /usr/lib/gcc/x86_64-w64-mingw32/*-posix/libgcc_s_seh-1.dll ~/winrb/bin
   cp /usr/lib/gcc/x86_64-w64-mingw32/*-posix/libstdc++-6.dll    ~/winrb/bin
   cp /usr/lib/gcc/x86_64-w64-mingw32/*-posix/libssp-0.dll       ~/winrb/bin
   cp ~/win_root/mingw64/bin/*.dll                               ~/winrb/bin
   ```

   This may be a bit fragile.  The paths for the first two DLLs are for Debian/Ubuntu, but could
   be somewhere else on other systems.

   Note that we need the `*-win32/` versions of `libgcc_s_seh-1.dll` and `libstdc++-6.dll` and not the
   `*-posix/` versions because we are using the mingw `libwinpthread-1.dll`.

   If you run `rb.exe` and it cannot find a DLL, it will tell you the first one that it cannot find.

Note that this same procedure is shown in the `release.yml` workflow.


## Troubleshooting:

* `rb: command not found`

    The problem is that you tried to run RevBayes but your computer doesn't know where the executable is. The easiest way is to add the directory in which you _installed_ RevBayes to your system path:

    ```
    export PATH=<prefix>/bin:$PATH  
    ```

    For example, if you set `-Dprefix=$HOME/Applications/revbayes`, you would do

    ```
    export PATH=$HOME/Applications/revbayes/bin:$PATH
    ```

* meson fails to configure

    You can look in `build/meson-logs/meson-log.txt` error log messages to help diagnose the problem.

* meson cannot find GTK2

    If GTK2 is not installed or not usable, this will produce an error message.

    Meson uses `pkg-config` to find most dependencies, and to find out what compiler and linker flags they need.

    If meson cannot find GTK2, you can:
    * run `pkg-config gtk+-2.0 --cflags` to see if pkg-config can find GTK2.
    * look in `build-gtk/meson-logs/meson-log.txt` for log messages to help diagnose the problem.

