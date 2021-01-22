#!/bin/bash

# Make a tarball of EPyGrAM tests

if [ "$1" == "-h" ]; then
    echo "Usage: mktar_tests.sh"
    exit
fi

VERSION=`grep __version__ epygram/__init__.py | awk '{print $3}' | awk -F "'" '{print $2}'`

tgz="$HOME/tmp/EPyGrAM-${VERSION}_tests.tgz"
tar -czf $tgz tests
ls -lh $tgz
cd ../..

