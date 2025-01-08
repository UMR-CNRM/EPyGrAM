#!/bin/bash

set -x

# FA
epygram cartoplot ../data/geometries/lambert_HS.fa -f SURFGEOPOTENTIEL
epygram cartoplot ../data/geometries/lambert_HS.fa -f SURFGEOPOTENTIEL -O /tmp/1.png
epygram cartoplot ../data/geometries/lambert_HS.fa -f SURFGEOPOTENTIEL --pm contourf --title contourf
epygram cartoplot ../data/geometries/lambert_HS.fa -f SURFGEOPOTENTIEL --pm contour --title contour
epygram cartoplot ../data/geometries/lambert_HS.fa -f SURFGEOPOTENTIEL --pm scatter --skw 's:20' --title scatter
epygram cartoplot ../data/geometries/lambert_HS.fa -f SURFGEOPOTENTIEL --cpyf RIVERS --title RIVERS
epygram cartoplot ../data/geometries/lambert_HS.fa -f SURFGEOPOTENTIEL --zoom 'lonmin=-80,lonmax=-70,latmin=-50,latmax=-45' --title zoom
epygram cartoplot ../data/geometries/lambert_HS.fa -f SURFGEOPOTENTIEL -m 0,1000 --title minmax
epygram cartoplot ../data/geometries/lambert_HS.fa -f SURFGEOPOTENTIEL -x '/,9.81' --title operation
