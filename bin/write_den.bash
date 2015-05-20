#!/bin/bash

# reduced version where the transformation happens in the python script

rm -v WORK/*iwfmt

cd WORK
echo -e "aoints\n 1\n" |$COLUMBUS/iwfmt.x > aoints.iwfmt 2> err.ignore

for file in *d1fl.* *d1trfl.*
do
    if [ -f $file ]; then
       echo $file
       echo -e "$file\n 1\n" |$COLUMBUS/iwfmt.x > $file.iwfmt 2> err.ignore
    fi
done
