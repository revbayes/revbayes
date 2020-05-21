#!/usr/bin/env bash
if [ -f "setup-x86_64.exe" ]; then
  ./setup-x86_64 -P `awk 'NR==1{printf \$1}{printf ",%s", \$1}' packagelist`
else
 echo "You must first copy cygwin 'setup-x86_64.exe' in the current folder."
fi
