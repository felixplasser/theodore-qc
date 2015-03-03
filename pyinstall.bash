#!/bin/bash

cd pyinstall || exit 2

for file in `ls ../bin/*.py` ../bin/theoinp
do
   pyinstaller -F $file
done
