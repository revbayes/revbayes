#!/bin/sh

set -e

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
    if [ -e ../../src/revlanguage/utils/GitVersion.cpp ] ; then
        cp ../../src/revlanguage/utils/GitVersion.cpp GitVersion_backup.cpp
    fi
    mv GitVersion.cpp ../../src/revlanguage/utils/

	./regenerate.sh $@
	cd build && cmake -G Ninja .
	ninja
	cd ..
	
    if [ -e  GitVersion_backup.cpp ] ; then
        cp GitVersion_backup.cpp ../../src/revlanguage/utils/GitVersion.cpp
        rm GitVersion_backup.cpp
    fi
fi
