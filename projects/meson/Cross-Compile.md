# Cross compiling with Meson

Cross-compiling means compiling executables on one system (e.g. Linux) that run on another system (e.g. Windows, Mac, ARM). Cross-compiling can be useful if you only have a Linux system, but want to create downloadable executables for another system.  One might also want to use Linux to compile ARM executables if the Linux system compiles faster.  The system the compiler runs on is called the "build" system, and the system that the executable will ultimately run on is called the "host" system.

## Linux to Windows

1. Download the cross-compiler and the build system
    ``` sh
    sudo apt-get install meson
    sudo apt-get install g++-mingw-w64 wine64       # Install cross-compiler and exe wrapper
    ```

2. Download RevBayes from our github repository.
    ``` sh
    git clone https://github.com/revbayes/revbayes.git revbayes
    ( cd revbayes/ ; git checkout development )     # Probably you want the development branch
    ```

3. Create a directory to hold windows libraries.
   ``` sh
   ./revbayes/projects/meson/make_winroot.sh       # Install windows boost in ~/win_root
   ```
   This script
   * Creates the "cross-file" `win64-cross.txt`
   * Create a pkg-config wrapper in ~/win_root/bin/pkg-config
   * Installs boost libraries in ~/win_root/mingw64/

4. Configure and compile 
   ``` sh
   ( cd revbayes/projects/meson/ ; ./generate.sh )
   meson build revbayes --prefix=$HOME/Applications/revbayes-w64 --cross-file=win64-cross.txt
   ninja -C build install
   ```

5. Copy the DLLs we use into the same directory as the revbayes executable
   ``` sh
   cp /usr/lib/gcc/x86_64-w64-mingw32/8.3-win32/libgcc_s_seh-1.dll $HOME/Applications/revbayes-64/bin
   cp /usr/lib/gcc/x86_64-w64-mingw32/8.3-win32/libstdc++-6.dll $HOME/Applications/revbayes-64/bin
   cp ~/win_root/mingw64/bin/*.dll $HOME/Applications/revbayes-64/bin
   ```
   You may have to change the path for libstdc++-6.dll and libgcc_s_seh-1.dll a bit, depending on
   your compiler version.

## Building RevStudio

For RevStudio we need to install win64 libraries for GTK2.  The libraries will end up in `~/win_root/mingw64/lib` with header files in `~/win_root/mingw64/include`.

   ``` sh
   cd ~/win_root

   for PKG in pango-1.42.4-2 \
              glib2-2.58.1-1 \
              cairo-1.16.0-1 \
              gdk-pixbuf2-2.38.0-2 \
              atk-2.30.0-1 \
              gtk2-2.24.32-3 \
              pcre-8.42-1 \
              fribidi-1.0.5-1 \
              pixman-0.36.0-1 \
              fontconfig-2.13.1-1 \
              freetype-2.9.1-1 \
              zlib-1.2.11-5 \
              bzip2-1.0.6-6 \
              libpng-1.6.36-1 \
              harfbuzz-2.2.0-1 \
              graphite2-1.3.12-1 \
              expat-2.2.6-1 \
              gettext-0.19.8.1-7 \
              libtiff-4.0.9-2 \
              jasper-2.0.14-1 \
              libjpeg-turbo-2.0.1-1 \
              xz-5.2.4-1 \
              libiconv-1.15-3 \
              libffi-3.2.1-3 \
              libthai-0.1.28-2 \
              libdatrie-0.2.12-1 \
              libwinpthread-git-7.0.0.5273.3e5acf5d-1 \
              ; do
      wget --no-verbose --show-progress http://repo.msys2.org/mingw/x86_64/mingw-w64-x86_64-${PKG}-any.pkg.tar.xz
      tar -Jxf mingw-w64-x86_64-${PKG}-any.pkg.tar.xz
   done
   ```

Then, you run meson with the `-Dstudio=true` flag:

``` sh
meson build-gtk revbayes -Dstudio=true -Dprefix=$HOME/Applications/revbayes-gui --cross-file=win64-cross.txt
ninja -C build-gtk install
```

It seems that statically linking the gtk libraries doesn't work, at least
when using the libraries from MSYS2.

## Installing Cross-compilers

The convention these days is that cross-compilers live in `/usr/bin/` and are prefixed with the description of the target architecture.  For example, `/usr/bin/x86_64-w64-mingw32-g++` is the compiler that targets 64-bit windows.

   ```
   apt-get install g++-mingw-w64
   apt-get install wine64
   ```

I think on homebrew, you can do:
   ```
   brew install mingw
   brew install wine
   ```
However, I have not tried this.


## Installing DLLs

   If your executables are not completely statically linked, you will have to
   include the DLLs that the executable needs.  If you put them in the same
   directory as the executable, then they will be found when you run it.

   It is possible to build rb.exe statically so that no DLLs are needed, but
   this doesn't seem possible for RevStudio.exe because static linking of
   GTK doesn't seem to be supported.  At the moment, static linking can be
   enabled or disabled with `-Dstatic=true` / `-Dstatic=false`.  This might
   be a bit of a hack.

   Libraries for GCC should be included with the compiler.  The exact location
   of these libraries depends on your linux/mac distribution.  On Debian, you
   might do something like this:

   ```
   cp /usr/lib/gcc/x86_64-w64-mingw32/8.3-win32/libgcc_s_seh-1.dll $OUTPUT/bin/
   cp /usr/lib/gcc/x86_64-w64-mingw32/8.3-win32/libstdc++-6.dll $OUTPUT/bin/
   cp ~/win_root/mingw64/bin/*.dll $OUTPUT/bin/
   ```

   It is possible that you might need slightly different libraries than the ones
   above, and its possible that you might need a 3rd one.

   For libraries related to BOOST and GTK, try running the executable under wine
   and see what DLLs it says are missing.  WINE claims DLLs are missing if they
   depend on another DLL that is missing, which is somewhat confusing.  But
   if you find all the DLLs that are missing and copy them to $OUTPUT/bin/ then
   the program should run under WINE.

   Currently RevStudio.exe and rb.exe both run under wine, but both stall after
   printing the welcome message and before the prompt.  However, they work
   in a windows VM, so this is probably a problem with WINE.



   