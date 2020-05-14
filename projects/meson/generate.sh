#!/bin/sh

set -e

SRC=../../src

./generate_help.sh

./generate_includes.sh

echo "
subdir('core')
subdir('revlanguage')
subdir('libs')
subdir('cmd')
" >> $SRC/meson.build

./generate_version_number.sh 
echo "Saving old GitVersion.cpp in projects/meson/GitVersion_backup.cpp"
if [ -e ../../src/revlanguage/utils/GitVersion.cpp ] ; then
    cp ../../src/revlanguage/utils/GitVersion.cpp GitVersion_backup.cpp
fi
echo "Copying current GitVersion.cpp to src/revlanguage/utils"
mv GitVersion.cpp ../../src/revlanguage/utils/

./generate_sources.sh core

./generate_sources.sh revlanguage

./generate_sources.sh libs

./generate_sources.sh cmd