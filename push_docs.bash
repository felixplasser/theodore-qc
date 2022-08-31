#!/bin/bash
# Udpate the documentation via the sourceforge server

cd doc || exit 1

# main html files
rsync -auv html/ fffelix@web.sourceforge.net:/home/project-web/theodore-qc/htdocs

# sphinx documentation
make html || exit 1
rsync -auv build/html/ fffelix@web.sourceforge.net:/home/project-web/theodore-qc/htdocs/docs
