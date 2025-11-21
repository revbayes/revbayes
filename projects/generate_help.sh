#!/bin/sh

PROJDIR=$(dirname $0)

SRC="$PROJDIR"/../src
HELP="$PROJDIR"/../help

echo "Updating RbHelpDatabase.cpp"

perl "$HELP"/md2help.pl "$HELP"/md/*.md > "$SRC"/core/help/RbHelpDatabase.cpp
