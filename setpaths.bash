#!/bin/bash
# source this file to set the paths

# Find the path of this script and set THEODIR automatically
ABSPATH=`readlink -e ${BASH_SOURCE[@]}`
export THEODIR=`dirname ${ABSPATH}`
# If this does not work set it manually and modify the next line
# export THEODIR=/user/plasserf/programs/TheoDORE/GIT

echo "THEODIR set to $THEODIR"

# set PATH and PYTHONPATH
export PATH=$THEODIR/bin:$PATH
export PYTHONPATH=$THEODIR

# add external packages here if they have not been installed in the default locations
#export PYTHONPATH=$HOME/programs/cclib/cclib-1.3.1/build/lib:$PYTHONPATH
