#!/bin/sh

set -e

PROJDIR=$(dirname $0)

SRC="$PROJDIR/../../src"

"$PROJDIR"/generate_includes.sh

echo "Defining variable 'help_pages' in 'src/meson.build'"
(
    cd "$SRC"
    echo "help_pages = files(["
    ls ../help/md/*.md | sed "s|^|'|" | sed "s|$|',|"
    echo "])"
) >> "$SRC"/meson.build

echo "
subdir('core')
subdir('revlanguage')
subdir('libs')
subdir('cmd')
subdir('help2yml')
" >> "$SRC"/meson.build

"$PROJDIR"/generate_version_number.sh
echo " Saving old GitVersion.cpp in projects/meson/GitVersion_backup.cpp"
if [ -e "$SRC"/revlanguage/utils/GitVersion.cpp ] ; then
    cp "$SRC"/revlanguage/utils/GitVersion.cpp GitVersion_backup.cpp
fi
echo " Copying current GitVersion.cpp to src/revlanguage/utils"
mv GitVersion.cpp "$SRC"/revlanguage/utils/

"$PROJDIR"/generate_sources.sh core

"$PROJDIR"/generate_sources.sh revlanguage

"$PROJDIR"/generate_sources.sh libs

"$PROJDIR"/generate_sources.sh cmd

"$PROJDIR"/generate_sources.sh help2yml
