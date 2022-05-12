#!/bin/bash
set -e

# Help
if [ "$1" == "-h" ]
then
    echo "You can pass the place to build as argument:"
    echo "- No argument will build in the directory of this script + /build/html"
    echo "- 'tmp' will build into a clone of https://github.com/UMR-CNRM/EPyGrAM-doc.git in a temporary directory in /tmp"
    echo "- or an existing directory (typically an EPyGrAM-doc repo)"
    exit
fi

# 0. Initialisations
WHERE_TO_BUILD=$1
BASEDIR=$(dirname $(readlink -f ${BASH_SOURCE[0]}))
export FOOTPRINT_DOCSTRINGS=2
DOC_FMT="html"
VERSION=`head -1 ../VERSION`
LAST_COMMIT=`git log -1 --oneline`

# 1. Output dir for the doc
if [ "$WHERE_TO_BUILD" == "" ]
then
    EPYGRAM_DOC=$BASEDIR/build/$DOC_FMT
    mkdir -p $EPYGRAM_DOC
elif [ "$WHERE_TO_BUILD" == "tmp" ]
then
    tmpdir=`mktemp -d -p /tmp`
    echo "Doc will be updated in '$tmpdir' as a temporary clone of https://github.com/UMR-CNRM/EPyGrAM-doc.git"
    cd $tmpdir
    git clone https://github.com/UMR-CNRM/EPyGrAM-doc.git
    EPYGRAM_DOC=$tmpdir/EPyGrAM-doc
    cd -
elif [ -d "$WHERE_TO_BUILD" ]
then
    EPYGRAM_DOC=$WHERE_TO_BUILD
    echo "Doc will be updated in existing '$EPYGRAM_DOC'"
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
  set +e
  bld=`which sphinx-build-3`
  if [ "$?" != "0" ]
  then
    bld=`which sphinx-build`
  fi
  set -e
  $bld -b $DOC_FMT . $EPYGRAM_DOC
cd ..

# 6. cleanings
find $EPYGRAM_DOC -name "*.ipynb" -print0 | xargs -0r rm
rm -rf $EPYGRAM_DOC/.doctrees

# Log
echo "==================================================="
echo "Doc built in $EPYGRAM_DOC for version: $LAST_COMMIT"

