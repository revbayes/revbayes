#!/bin/sh

SRC=../../src

./generate_includes.sh

echo "
subdir('core')
subdir('revlanguage')
subdir('libs')
subdir('cmd')
" >> $SRC/meson.build

./generate_sources.sh core

./generate_sources.sh revlanguage

./generate_sources.sh libs

./generate_sources.sh cmd

./generate_version_number.sh 
echo "Saving old GitVersion.cpp in projects/meson/GitVersion_backup.cpp"
cp ../../src/revlanguage/utils/GitVersion.cpp GitVersion_backup.cpp
echo "Copying current GitVersion.cpp to src/revlanguage/utils"
mv GitVersion.cpp ../../src/revlanguage/utils/
