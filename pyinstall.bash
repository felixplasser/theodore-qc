#!/bin/bash
# call this on vanadium

cd pyinstall || exit 2

for file in `ls ../bin/*.py` ../bin/theoinp
do
#   pyinstaller -F $file
   pyinstaller --strip -F $file
done

# strip could make smaller files. But does this work?
#cd dist
#strip *
