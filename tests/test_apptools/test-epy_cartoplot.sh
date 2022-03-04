#!/bin/bash

python=`which python3`

if [ "$?" != "0" ]
then
  echo "ERROR: do not know which 'python3' to use"
  exit 1
fi

set -x

# FA
$python ../../apptools/epy_cartoplot.py ../data/geometries/lambert_HS.fa -f SURFGEOPOTENTIEL
$python ../../apptools/epy_cartoplot.py ../data/geometries/lambert_HS.fa -f SURFGEOPOTENTIEL -O /tmp/1.png
$python ../../apptools/epy_cartoplot.py ../data/geometries/lambert_HS.fa -f SURFGEOPOTENTIEL --pm contourf --title contourf
$python ../../apptools/epy_cartoplot.py ../data/geometries/lambert_HS.fa -f SURFGEOPOTENTIEL --pm contour --title contour
$python ../../apptools/epy_cartoplot.py ../data/geometries/lambert_HS.fa -f SURFGEOPOTENTIEL --pm scatter --skw 's:20' --title scatter
$python ../../apptools/epy_cartoplot.py ../data/geometries/lambert_HS.fa -f SURFGEOPOTENTIEL --cpyf RIVERS --title RIVERS
$python ../../apptools/epy_cartoplot.py ../data/geometries/lambert_HS.fa -f SURFGEOPOTENTIEL --zoom 'lonmin=-80,lonmax=-70,latmin=-50,latmax=-45' --title zoom
$python ../../apptools/epy_cartoplot.py ../data/geometries/lambert_HS.fa -f SURFGEOPOTENTIEL -m 0,1000 --title minmax
$python ../../apptools/epy_cartoplot.py ../data/geometries/lambert_HS.fa -f SURFGEOPOTENTIEL -x '/,9.81' --title operation
