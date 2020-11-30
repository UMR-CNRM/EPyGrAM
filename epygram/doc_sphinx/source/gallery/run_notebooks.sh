#!/usr/bin/bash

arg=$1

if [ "$arg" == "in_place" ]
then
  cmdline="--inplace"
elif [ "$arg" == "check_errors" ]
then
  cmdline='--stdout'
else
  echo "Arguments:"
  echo " in_place : run the notebooks and save the status in place"
  echo " check_errors : just check that no error arise"
  exit
fi

for f in `ls */*.ipynb`
do
  jupyter nbconvert --execute $cmdline $f >> /dev/null
  exit
done
