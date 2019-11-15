#!/bin/bash
export FOOTPRINT_DOCSTRINGS=2
fmt="html"
if [ "$1" != "" ]; then
  fmt=$1
fi
echo "> Create cmaps..."
dynamic/make_cmaps_png_for_doc.py
echo "> List external dependancies..."
../../list_external_dependancies.py
cd source
echo "> Build cheatsheet..."
pdflatex cheatsheet.tex
rm -f cheatsheet.log cheatsheet.aux
cp -f cheatsheet.pdf _static/.
mv -f cheatsheet.pdf ../html/_downloads/cheatsheet.pdf
echo "> Build Sphinx doc..."
sphinx-build -b $fmt . ../$fmt
cd ..
if [ "$fmt" == "latex" ]; then
  cd $fmt
  make
fi
