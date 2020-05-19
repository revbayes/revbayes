Building RevBayes with cmake
############################

You build RevBayes by running a build script appropriate for your
operating system:

Script              Operating system
==================  ==============================================
build.sh            Unix, Linux and Mac OS X

Example for a Linux system:
$ ./build.sh

The build script produces the RevBayes console executable 'rb'

Additional Options
===================  =============================================
-boost <true|false>  (re)compile Boost libs. Default=true
-mac   <true|false>  build for OS X 10.6 and higher. Default=false
-win   <true|false>  build for a Windows system. Default=false
-mpi   <true|false>  build the MPI version 'rb-mpi'. Default=false

More Information
================
The build script invokes the regenerate script to build the boost
libraries and preparing things for cmake. Then it runs cmake, and
finally it uses make to compile a set of executables.

To avoid using the parallel version of make, you will need to
modify the build.sh script by deleting the -j 4 flag:
$ make -j 4
