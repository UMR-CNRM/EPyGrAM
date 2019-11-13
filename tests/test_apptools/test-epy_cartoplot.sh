#!/bin/bash
set -x

# FA
epy_cartoplot.py ../data/geometries/lambert_HS.fa -f SURFGEOPOTENTIEL
epy_cartoplot.py ../data/geometries/lambert_HS.fa -f SURFGEOPOTENTIEL -O /tmp/1.png
epy_cartoplot.py ../data/geometries/lambert_HS.fa -f SURFGEOPOTENTIEL --pm contourf --title contourf
epy_cartoplot.py ../data/geometries/lambert_HS.fa -f SURFGEOPOTENTIEL --pm contour --title contour
epy_cartoplot.py ../data/geometries/lambert_HS.fa -f SURFGEOPOTENTIEL --pm scatter --skw 'pointsize:20' --title scatter
epy_cartoplot.py ../data/geometries/lambert_HS.fa -f SURFGEOPOTENTIEL --cpyf RIVERS --title RIVERS
epy_cartoplot.py ../data/geometries/lambert_HS.fa -f SURFGEOPOTENTIEL --zoom 'lonmin=-80,lonmax=-70,latmin=-50,latmax=-45' --title zoom
epy_cartoplot.py ../data/geometries/lambert_HS.fa -f SURFGEOPOTENTIEL -m 0,1000 --title minmax
epy_cartoplot.py ../data/geometries/lambert_HS.fa -f SURFGEOPOTENTIEL -x '/,9.81' --title operation