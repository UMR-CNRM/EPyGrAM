#!/bin/bash

WHERE_TO_BUILD=$1

BASEDIR=$(dirname $(readlink -f ${BASH_SOURCE[0]}))

export FOOTPRINT_DOCSTRINGS=2
DOC_FMT="html"
VERSION=`head -1 ../../VERSION`

# 1. Output dir for the doc
if [ "$WHERE_TO_BUILD" == "" ]
then
    EPYGRAM_DOC=$BASEDIR/build/$DOC_FMT
    mkdir -p $EPYGRAM_DOC
elif [ "WHERE_TO_BUILD" == "tmp" ]
then
    tmpdir=`mktemp -d -p /tmp`
    if [ "$?" != "0" ]
    then
        echo "Error creating tmpdir in /tmp"
        exit 1
    else
        echo "Build in tmpdir=$tmpdir"
    fi
    echo "Doc will be updated in '$tmpdir' as a temporary clone of https://github.com/UMR-CNRM/EPyGrAM-doc.git"
    cd $tmpdir
    git clone https://github.com/UMR-CNRM/EPyGrAM-doc.git
    EPYGRAM_DOC=$tmpdir/EPyGrAM-doc
    cd -
else
    echo "Error: unrecognized argument #1: '$WHERE_TO_BUILD', must be 'tmp' if argument is provided"
    echo "Exit"
    exit 1
fi

# 2. Cmaps
echo "> Create cmaps..."
python3 dynamic/make_cmaps_png_for_doc.py -d $EPYGRAM_DOC

# 3. Dependancies
echo "> List external dependancies..."
python3 dynamic/list_external_dependancies.py

# Source compilation
cd source

  # 4. Cheatsheet
  echo "> Build cheatsheet..."
  pdflatex cheatsheet.tex
  if [ ! -d _static ]
  then
    mkdir _static
  fi
  cp -f cheatsheet.pdf _static/.

  # 5. Sphinx
  echo "> Build Sphinx doc..."
  bld=`which sphinx-build-3`
  if [ "$?" != "0" ]
  then
    bld=`which sphinx-build`
  fi
  $bld -b $DOC_FMT . $EPYGRAM_DOC
cd ..

# 6. cleanings
find $EPYGRAM_DOC -name "*.ipynb" -print0 | xargs -0r rm
rm -rf $EPYGRAM_DOC/.doctrees

# Log
echo "=============================================="
echo "Doc built for version $VERSION in $EPYGRAM_DOC"

