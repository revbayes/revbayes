#!/bin/sh

# see directions in /projects/RevBayes_Win_Installation_Instruction.txt

if [ "$1" = "clean" ]
then
	rm -rf build
else
if [ ! -d build ]; then
	mkdir build
fi

    #################
    # generate git version number
    ./generate_version_number.sh
    cp ../../src/revlanguage/utils/GitVersion.cpp GitVersion_backup.cpp
    mv GitVersion.cpp ../../src/revlanguage/utils/

	./regenerate.sh -win true -boost false $@
	cd build 
	cmake -v -DCMAKE_TOOLCHAIN_FILE=../mingw64_toolchain.cmake .
	make -j 4
	cd ..
	
    cp GitVersion_backup.cpp ../../src/revlanguage/utils/GitVersion.cpp
    rm GitVersion_backup.cpp
fi
