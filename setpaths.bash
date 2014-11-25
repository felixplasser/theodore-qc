#!/bin/bash
# source this file to set the paths

export THEODIR=/export/home/fplasser/programs/TheoDORE/GIT
export PATH=$THEODIR/bin:$PATH
export PYTHONPATH=$THEODIR/lib

# if the external cclib package is used:
export PYTHONPATH=/export/home/fplasser/programs/cclib/cclib-1.3/build/lib.linux-x86_64-2.7:$PYTHONPATH
