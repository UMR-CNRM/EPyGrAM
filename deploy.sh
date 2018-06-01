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
else
    version='-'$version
fi
mkdoc=$2
if [ "$mkdoc" == "mkdoc" ]; then
    cd epygram/doc_sphinx
    ./mk_html_doc.sh
    cd ../..
fi


# Filter
to_exclude=''
for elem in playground tests versioning.txt deploy.sh apptools/*.pyc site/arpifs4py/libs4py_*.so epygram/doc_sphinx/source/gallery/inputs *__pycache__*
do
  to_exclude="$to_exclude --exclude $elem"
done
echo $to_exclude
no_pyc='--exclude *.pyc'


# Make temp dir, copy, go therein
tmp=`mktemp`
if [ "$tmp" == "" ]; then
    echo "Problem in tmp dir generation. Exit."
    exit
fi
rm -f $tmp; mkdir $tmp
rsync -avL * $tmp/ $to_exclude
here=`pwd`
cd $tmp
echo "Tempdir:" `pwd`


# Rsync
LOC_SXCOOPE="sxcoope1:~mary/sync_epygram/EPyGrAM$version/"
rsync -av * $LOC_SXCOOPE

#LOC_VXDEV64="vxdev64:~mary/EPyGrAM$version/" # vxdev64: development server @ CNRM (OS updates)
#rsync -av * $LOC_VXDEV64 $no_pyc

rm site/arpifs4py/libs4py.so

LOC_PAGRE="pagre:~mary/public/EPyGrAM$version/"
rsync -av * $LOC_PAGRE $no_pyc

rm -rf site/epyweb

LOC_BFX="beaufix:~mary/public/EPyGrAM$version/"
rsync -av * $LOC_BFX $no_pyc
LOC_PLX="prolix:~mary/public/EPyGrAM$version/"
rsync -av * $LOC_PLX $no_pyc
echo ""
echo "==> deployment done on: $LOC_SXCOOPE $LOC_VXDEV64 $LOC_PAGRE $LOC_BFX $LOC_PLX"
echo "libs4py.so to be linked on beaufix/prolix (~mary/deploy_epygram_finalize.sh) and pagre"

# Come back and clean
cd $here
rm -rf $tmp

