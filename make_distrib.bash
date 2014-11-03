#!/bin/bash
# Create a new distribution version

echo "Syntax: make_distrib.bash [Version]"

if [ $# -eq 0 ]
then
   echo "Please enter version number!"
fi
