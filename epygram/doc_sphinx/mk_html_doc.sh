#!/bin/bash
export FOOTPRINT_DOCSTRINGS=2
fmt="html"
if [ "$1" != "" ]; then
  fmt=$1
fi
cd source
sphinx-build -b $fmt . ../$fmt
cd ..
if [ "$fmt" == "latex" ]; then
  cd $fmt
  make
fi
