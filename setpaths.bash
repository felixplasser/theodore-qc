#!/bin/bash
# source this file to set the paths

export THEODIR=$HOME/programs/TheoDORE/GIT
export PATH=$THEODIR/bin:$PATH
export PYTHONPATH=$THEODIR/lib

# add external packages here if they have not been installed in the default locations
export PYTHONPATH=$HOME/programs/cclib/cclib-1.3/build/lib.linux-x86_64-2.7:$PYTHONPATH
