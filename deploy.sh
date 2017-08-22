#!/bin/bash

# Script de deploiement d'EPyGrAM sur /home/common/epygram et BULL
# Reste manuel: déploiement .tar sur redmine et hendrix (à automatiser)

# Parse args
if [ "$1" == "-h" ]; then
    echo "Usage: deploy.sh version [mkdoc]"
    echo "version being e.g. 'dev' or '1.1.8'"
    echo "if mkdoc is present, build doc before deployment"
    exit
fi
version=$1
if [ "$version" == "" ]; then
    echo "Need to provide version as argument !"
    exit
fi
mkdoc=$2
if [ "$mkdoc" == "mkdoc" ]; then
    cd epygram/doc_sphinx
    ./mk_html_doc.sh
    cd ../..
fi

# Make temp dir, copy, go therein
tmp=`mktemp`
if [ "$tmp" == "" ]; then
    echo "Problem in tmp dir generation. Exit."
    exit
fi
rm -f $tmp; mkdir $tmp
cp -rL * $tmp/.
here=`pwd`
cd $tmp
echo "Tempdir:" `pwd`

# Filter
rm -rf playground
rm -rf tests
rm -f versioning.txt
rm -f deploy.sh
rm -f apptools/*.pyc
rm site/arpifs4py/libs4py_*.so

# Rsync
rsync -av * sxcoope1:~mary/sync_epygram/EPyGrAM.$version/
rm site/arpifs4py/libs4py.so
rm -rf site/epyweb
rsync -av * beaufix:~mary/public/EPyGrAM.$version/
rsync -av * prolix:~mary/public/EPyGrAM.$version/
echo "libs4py.so to be linked on beaufix/prolix (~mary/deploy_epygram_finalize.sh)"

# Come back and clean
cd $here
rm -rf $tmp

