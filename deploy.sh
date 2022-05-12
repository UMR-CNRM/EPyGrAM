#!/bin/bash

# Script de deploiement d'EPyGrAM sur /home/common/epygram (CNRM) et autres machines

# Parse args
if [ "$1" == "-h" ]; then
    echo "Usage: deploy.sh [VERSION]"
    echo "<VERSION> being the distant label, e.g. 'dev'"
    echo "If no VERSION is provided, the numbered version found in epygram/__init__.py is used."
	echo "The distant installation is labelled EPyGrAM-<VERSION>"
    exit
fi
VERSION=$1
if [ "$VERSION" == "" ]; then
    VERSION=`head -1 VERSION`
fi
EPYGRAM_DIR="public/EPyGrAM/$VERSION"


# Platforms to push onto
sxcoope1=1  # from which are synchronised all CNRM workstations
vxdev64=0  # vxdev64: development server @ CNRM (OS updates)
sotrtm33sidev=1  # COMPAS server, from which it is replicated onto the others
belenos=1
taranis=1


# Filters
to_exclude4all=''
for elem in playground tests deploy.sh mktar.sh apptools/*.pyc site/arpifs4py/libs4py_*.so docs/source/gallery/inputs *.ipynb_checkpoints* docs/source/gallery/inputs/* gallery/inputs/*
do
  to_exclude4all="$to_exclude4all --exclude $elem"
done
no_source_doc=''
for elem in docs/source/* docs/build/* docs/Makefile docs/mk_html_doc.sh
do
  no_source_doc="$no_source_doc --exclude $elem"
done
no_pyc='--exclude *.pyc --exclude *__pycache__*'
no_doc='--exclude docs/*'
no_libs4py='--exclude site/arpifs4py/libs4py.so'
no_arpifs4py='--exclude site/arpifs4py'
no_epyweb='--exclude site/epyweb'

# Filters specific to platforms
to_exclude4sxcoope1="$to_exclude4all"
to_exclude4vxdev64="$to_exclude4all $no_pyc $no_source_doc"
to_exclude4sidev="$to_exclude4all $no_pyc $no_arpifs4py"
to_exclude4bull="$to_exclude4all $no_pyc $no_source_doc $no_libs4py $no_epyweb"


# Rsync
logger="EPyGrAM:$VERSION deployed on:\n"
echo "------------------------------------------------------"
if [ "$vxdev64" == 1 ]; then
  echo "...vxdev64..."
  rsync -avL * vxdev64:$EPYGRAM_DIR $to_exclude4vxdev64
  logger="$logger - vxdev64\n"
fi
echo "------------------------------------------------------"
if [ "$sxcoope1" == 1 ]; then
  echo "...sxcoope1..."
  rsync -avL * sxcoope1:$EPYGRAM_DIR $to_exclude4sxcoope1
  logger="$logger - sxcoope1\n"
fi
echo "------------------------------------------------------"
if [ "$sotrtm33sidev" == 1 ]; then
  echo "...sotrtm33-sidev..."
  rsync -avL * sotrtm33-sidev:$EPYGRAM_DIR $to_exclude4sidev
  logger="$logger - sotrtm33-sidev\n"
fi
echo "------------------------------------------------------"
if [ "$belenos" == 1 ]; then
  echo "...belenos..."
  rsync -avL * belenos:$EPYGRAM_DIR $to_exclude4bull
  logger="$logger - belenos\n"
fi
echo "------------------------------------------------------"
if [ "$taranis" == 1 ]; then
  echo "...taranis..."
  rsync -avL * taranis:$EPYGRAM_DIR $to_exclude4bull
  logger="$logger - taranis\n"
fi


# Log final
echo "------------------------------------------------------"
echo -e $logger
echo "Don't forget to link *libs4py.so* on necessary machines (supercomputers) !"

