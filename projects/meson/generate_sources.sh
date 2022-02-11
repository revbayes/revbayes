#!/bin/sh

DIR=$1 

echo "Defining variable '${DIR}_sources' in src/${DIR}/meson.build"

SRC=../../src

cd $SRC/$DIR

echo "${DIR}_sources = files([" > meson.build
find . -name '*.cpp' |
    sed "s|^|'|" |
    sed "s|$|',|" |
    grep -v '^\.' |
    grep -v 'main.cpp' >> meson.build
find . -name '*.c' |
    sed "s|^|'|" |
    sed "s|$|',|" |
    grep -v '^\.' >> meson.build
echo "])" >> meson.build
echo >> meson.build
