#!/bin/bash

# Make a tarball of EPyGrAM doc

if [ "$1" == "-h" ]; then
    echo "Usage: mktar_doc.sh"
    exit
fi

VERSION=`head -1 VERSION`

cd docs
tgz="$HOME/tmp/EPyGrAM-${VERSION}_doc.tgz"
tar -czf $tgz build
ls -lh $tgz
cd ..

