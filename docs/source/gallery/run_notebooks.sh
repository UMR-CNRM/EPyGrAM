#!/usr/bin/bash

arg=$1

if [ "$arg" == "in_place" ]
then
  args="--execute --inplace"
elif [ "$arg" == "check_errors" ]
then
  args='--execute --stdout'
elif [ "$arg" == "clear_output" ]
then
  args='--clear-output --inplace --ClearOutputPreprocessor.enabled=True'
else
  echo "Arguments:"
  echo " in_place : run the notebooks and save the status in place"
  echo " check_errors : just check that no error arise"
  echo " clear_output : clear outputs to make the notebooks way lighter (to be done before committing)"
  exit
fi

for f in `ls */*.ipynb`
do
  jupyter nbconvert $args $f >> /dev/null
done
