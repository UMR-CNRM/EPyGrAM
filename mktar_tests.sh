#!/bin/bash

# Make a tarball of EPyGrAM tests

if [ "$1" == "-h" ]; then
    echo "Usage: mktar_tests.sh"
    exit
fi

VERSION=`head -1 VERSION`

tgz="$HOME/tmp/EPyGrAM-${VERSION}_tests.tgz"
tar -czf $tgz tests
ls -lh $tgz

