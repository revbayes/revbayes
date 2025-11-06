#!/bin/sh

PROJDIR=$(dirname $0)

SRC="$PROJDIR"/../../src

cd "$SRC"

echo "Defining variable 'src_inc' in 'src/meson.build'"

echo "src_inc = include_directories([" > meson.build
find libs revlanguage core cmd help2yml -name '*.h*' |
    sed "s|/[^/]*$||g" |
    sed "s|libs/\([^/]*\)/.*|libs/\1|" |
    sort |
    uniq |
    sed "s|^|'|" |
    sed "s|$|',|" >> meson.build
echo "])" >> meson.build
echo >> meson.build
