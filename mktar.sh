#!/bin/bash

# Make a tarball of EPyGrAM (without tests nor doc, to be tared separately)

if [ "$1" == "-h" ]; then
    echo "Usage: mktar.sh [/path/to/vortex]"
    exit
fi

VORTEX="/home/common/sync/vortex/vortex"
if [ "$1" != "-h" ] && [ "$1" != "" ]; then
    VORTEX=$1
fi

VERSION=`head -1 VERSION`

# Filters
to_exclude=''
for elem in tests deploy.sh mktar.sh apptools/*.pyc site/arpifs4py/libs4py_*.so docs/build docs/source/gallery/inputs docs/source/gallery/outputs *__pycache__* */.ipynb_checkpoints/* site/usevortex.py
do
  to_exclude="$to_exclude --exclude $elem"
done
echo $to_exclude
no_pyc='--exclude *.pyc'
no_pycache='--exclude *__pycache__*'
no_libs4py='--exclude site/arpifs4py/libs4py.so'

# Make temp dir, copy, go therein
tmp=`mktemp`
if [ "$tmp" == "" ]; then
    echo "Problem in tmp dir generation. Exit."
    exit
fi
rm -f $tmp; mkdir $tmp
rsync -avL * $tmp/ $to_exclude $no_pyc $no_pycache $no_libs4py
here=`pwd`
cd $tmp
echo "Tempdir:" `pwd`

# Add necessary vortex site packages
vortex_dependancies=${VORTEX}/site
for package in "bronx" "footprints" "taylorism"
do
  rsync -av $no_pyc $no_pycache $vortex_dependancies/$package site/
done

# Tar
tgz="$HOME/tmp/EPyGrAM-$VERSION.tgz"
tar -czf $tgz *
ls -lh $tgz

# Come back and clean
cd $here
rm -rf $tmp

