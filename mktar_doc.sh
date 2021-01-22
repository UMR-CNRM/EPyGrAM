#!/bin/bash

# Make a tarball of EPyGrAM doc

if [ "$1" == "-h" ]; then
    echo "Usage: mktar_doc.sh"
    exit
fi

VERSION=`grep __version__ epygram/__init__.py | awk '{print $3}' | awk -F "'" '{print $2}'`

cd epygram/doc_sphinx
tgz="$HOME/tmp/EPyGrAM-${VERSION}_doc.tgz"
tar -czf $tgz html
ls -lh $tgz
cd ../..

