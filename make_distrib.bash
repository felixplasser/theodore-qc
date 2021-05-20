#!/bin/bash
# Create a new distribution version

# TODO before:
#   - Check version number in theo_header.py
#   - git tag set by this script

echo "Syntax: make_distrib.bash [Version]"
echo

if [ $# -eq 0 ]
then
   echo "Please enter version number!"
   exit 1
fi

SDIR=`pwd`

LTDIR="TheoDORE_$1"
TDIR="$SDIR/../Versions/$LTDIR"
echo "Creating new version in $TDIR"

git tag -a v$1

mkdir $TDIR || exit 2

cp README COPYRIGHT.txt LICENSE.txt $TDIR

sed "s/GIT/$LTDIR/" setpaths.bash > $TDIR/setpaths.bash
sed "s/GIT/$LTDIR/" setpaths.csh > $TDIR/setpaths.csh

cp -r bin $TDIR
cp -r EXAMPLES $TDIR

mkdir $TDIR/theodore
cp theodore/*.py $TDIR/theodore

# cclib as used by TheoDORE
cp -r external/cclib/cclib $TDIR
cp external/cclib/LICENSE $TDIR/cclib
echo "Removing cclib binary files"
rm -r $TDIR/cclib/__pycache__ $TDIR/cclib/*/__pycache__

# remove extra files in EXAMPLE directory
rm -r $TDIR/EXAMPLES/*/RUN
rm -r $TDIR/EXAMPLES/*/__pycache__

# create tar with shorter relative paths
cd $SDIR/../Versions
tar -czf $LTDIR.tar.gz $LTDIR
