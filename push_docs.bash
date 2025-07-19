#!/bin/bash
# Udpate the documentation via the sourceforge server

cd doc || exit 1

# main html files
rsync -auv html/ fffelix@web.sourceforge.net:/home/project-web/theodore-qc/htdocs

# sphinx documentation
../bin/print_doc.py > source/keywords.rst || exit 1
make html || exit 1
rsync -auv build/html/ fffelix@web.sourceforge.net:/home/project-web/theodore-qc/htdocs/docs

# push from docs2
echo
echo " *** Now run: "
echo " git checkout docs2 "
echo " ./push_github_docs.bash "
