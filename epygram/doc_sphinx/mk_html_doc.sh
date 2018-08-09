#!/bin/bash
export FOOTPRINT_DOCSTRINGS=2
fmt="html"
if [ "$1" != "" ]; then
  fmt=$1
fi
dynamic/make_cmaps_png_for_doc.py
../../list_external_dependancies.py
cd source
sphinx-build -b $fmt . ../$fmt
cd ..
if [ "$fmt" == "latex" ]; then
  cd $fmt
  make
fi
