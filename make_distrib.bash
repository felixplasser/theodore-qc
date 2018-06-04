#!/bin/bash
# Create a new distribution version

# TODO before:
#   - Check version number in theo_header.py

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

if [ -d $TDIR ]
then
    echo " $TDIR already exists!"
    exit 2
fi

mkdir $TDIR

cp README COPYRIGHT.txt LICENSE.txt $TDIR

sed "s/GIT/$LTDIR/" setpaths.bash > $TDIR/setpaths.bash
sed "s/GIT/$LTDIR/" setpaths.csh > $TDIR/setpaths.csh

cp -r bin $TDIR
cp -r EXAMPLES $TDIR

mkdir $TDIR/lib
cp lib/*.py $TDIR/lib

# cclib as used by TheoDORE
echo "Removing cclib binary files"
rm -v external/cclib/cclib/*.pyc external/cclib/cclib/*/*.pyc
cp -r external/cclib/cclib $TDIR/lib
# also copy the source because of LGPL
mkdir $TDIR/external
cd $SDIR/external
tar -czf $TDIR/external/cclib.tgz cclib/ANNOUNCE cclib/CHANGELOG cclib/INSTALL cclib/LICENSE cclib/logo_for_ccl.svg cclib/logo.png cclib/manifest.py cclib/README.md cclib/setup.py cclib/src/ cclib/test/ cclib/THANKS

# create tar with shorter relative paths
cd $SDIR/../Versions
tar -czf $LTDIR.tgz $LTDIR
