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

LTDIR="TheoDORE_$1"
TDIR="../$LTDIR"
echo "Creating new version in $TDIR"

if [ -d $TDIR ]
then
    echo " $TDIR already exists!"
    exit 2
fi

mkdir $TDIR

cp COPYRIGHT.txt LICENSE.txt $TDIR

sed "s/GIT/$LTDIR/" setpaths.bash > $TDIR/setpaths.bash
sed "s/GIT/$LTDIR/" setpaths.csh > $TDIR/setpaths.csh

cp -r bin $TDIR
cp -r EXAMPLES $TDIR

mkdir $TDIR/lib
cp lib/*.py $TDIR/lib

cd ..
tar -czf $LTDIR.tgz $LTDIR
