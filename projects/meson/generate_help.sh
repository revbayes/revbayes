#!/bin/sh


SRC=../../src
HELP=../../help

echo "Updating RbHelpDatabase.cpp"

perl $HELP/md2help.pl $HELP/md/*.md > $SRC/core/help/RbHelpDatabase.cpp
