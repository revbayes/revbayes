#!/bin/sh

echo "Generating version numbers in GitVersion.cpp"
echo "#include \"GitVersion.h\"" > GitVersion.cpp
#echo "const char *build_git_sha = \"$(git rev-parse HEAD)\";" >> GitVersion.cpp
echo "const char *build_git_sha = \"$(git describe --abbrev=6 --always)\";" >> GitVersion.cpp
echo "const char *build_date = \"$(date)\";" >> GitVersion.cpp
echo "const char *build_git_branch = \"$(git name-rev --name-only HEAD)\";" >> GitVersion.cpp
