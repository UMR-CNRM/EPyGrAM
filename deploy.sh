#!/bin/bash

# Script de deploiement d'EPyGrAM sur /home/common/epygram et BULL

# Parse args
if [ "$1" == "-h" ]; then
    echo "Usage: deploy.sh version [mkdoc]"
    echo "version being the distant label, e.g. 'dev' or '1.1.8'"
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
sxcoope1=1
vxdev64=0
pagre=1
bullx=1


# Filter
to_exclude4all=''
for elem in playground tests versioning.txt deploy.sh mktar.sh apptools/*.pyc site/arpifs4py/libs4py_*.so epygram/doc_sphinx/source/gallery/inputs *.ipynb_checkpoints*
do
  to_exclude4all="$to_exclude4all --exclude $elem"
done
no_source_doc=''
for elem in epygram/doc_sphinx/source/* epygram/doc_sphinx/build/* epygram/doc_sphinx/Makefile epygram/doc_sphinx/mk_html_doc.sh
do
  no_source_doc="$no_source_doc --exclude $elem"
done
no_pyc='--exclude *.pyc --exclude *__pycache__*'
no_doc='--exclude epygram/doc_sphinx/*'
no_libs4py='--exclude site/arpifs4py/libs4py.so'
no_arpifs4py='--exclude site/arpifs4py'
no_epyweb='--exclude site/epyweb'


# Rsync
# sxcoope1
if [ "$sxcoope1" == 1 ]; then
  echo "------------------------------------------------------"
  echo "...sxcoope1..."
  LOC_SXCOOPE="sxcoope1:~mary/sync_epygram/EPyGrAM$version/"  # from which are synchronised all CNRM workstations
  rsync -avL * $LOC_SXCOOPE $to_exclude4all
fi
# vxdev64
if [ "$vxdev64" == 1 ]; then
  echo "------------------------------------------------------"
  echo "...vxdev64..."
  LOC_VXDEV64="vxdev64:~mary/EPyGrAM$version/"  # vxdev64: development server @ CNRM (OS updates)
  rsync -avL * $LOC_VXDEV64 $to_exclude4all $no_pyc $no_source_doc
fi
# pagre
if [ "$pagre" == 1 ]; then
  echo "------------------------------------------------------"
  echo "...pagre..."
  LOC_PAGRE="pagre:~mary/public/EPyGrAM$version/"  # COMPAS server
  rsync -avL * $LOC_PAGRE $to_exclude4all $no_pyc $no_arpifs4py
fi
# bullx
if [ "$bullx" == 1 ]; then
  echo "------------------------------------------------------"
  echo "...bullx..."
  bull_exclude="$to_exclude4all $no_pyc $no_source_doc $no_libs4py $no_epyweb"
  bull_public="~mary/public/EPyGrAM$version/"
  rsync -avL * "beaufix:$bull_public" $bull_exclude
  rsync -avL * "prolix:$bull_public" $bull_exclude
fi


# Summary
echo "------------------------------------------------------"
if [ "$sxcoope1" == 1 ]; then
  echo "=> deployed on sxcoope1"
fi
if [ "$vxdev64" == 1 ]; then
  echo "=> deployed on vxdev64"
fi
if [ "$pagre" == 1 ]; then
  echo "=> deployed on pagre"
  echo "!!! Deactivate arpifs4py formats there !!!"
fi
if [ "$bullx" == 1 ]; then
  echo "=> deployed on beaufix & prolix"
  echo "   libs4py.so to be linked there (~mary/deploy_epygram_finalize.sh)"
fi

