#!/bin/bash
# Create a new distribution version

# TODO before:
#   - Check version number in theo_header.py and doc/source/conf.py
#   - doc/README up to date?
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

mkdir $TDIR || exit 2

git tag -a v$1

cp README.rst COPYRIGHT.txt LICENSE.txt $TDIR

sed "s/GIT/$LTDIR/" setpaths.bash > $TDIR/setpaths.bash
sed "s/GIT/$LTDIR/" setpaths.csh > $TDIR/setpaths.csh

cp -r bin $TDIR
mkdir $TDIR/EXAMPLES
cp -r ../EXAMPLES/STANDARD $TDIR/EXAMPLES
cp -r ../EXAMPLES/CCLIB    $TDIR/EXAMPLES
cp -r ../EXAMPLES/UTILS    $TDIR/EXAMPLES
cp -r ../EXAMPLES/EXTRA    $TDIR/EXAMPLES
rm -r $TDIR/EXAMPLES/*/*/RUN

cp -r theodore $TDIR

# cclib as used by TheoDORE
cp -r external/cclib/cclib $TDIR
cp external/cclib/LICENSE $TDIR/cclib

# periodictable
cp -r external/periodictable/periodictable $TDIR
cp external/periodictable/LICENSE.txt $TDIR/periodictable

# colt as used by TheoDORE
cp -r external/colt/colt $TDIR
cp external/colt/LICENSE $TDIR/colt

echo "Removing binary pyc files"
find $TDIR -name '*pyc' -exec rm {} \;

# create tar with shorter relative paths
cd $SDIR/../Versions
tar -czf $LTDIR.tar.gz $LTDIR
