#!/bin/bash
# Preparation of documentation:

# Change relevant files in master branch

git pull origin main || exit 1

cd doc
make html || exit 1
cd ..
/bin/cp -r doc/html/* docs
/bin/cp -r doc/build/html/* docs/docs

# submit docs to GIT
git commit -m 'doc update' docs
git push origin docs2
