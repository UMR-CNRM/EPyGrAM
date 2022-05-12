#!/usr/bin/bash

HOST=ftp.umr-cnrm.fr
PORT=21
LOGIN=epygram  # read-only anonymous user
PASSWORD=EpygraM  # with anonymous password
INPUTS_DIR=inputs

if [ ! -d "$INPUTS_DIR" ]
then
    mkdir $INPUTS_DIR
fi

cd $INPUTS_DIR
ftp -i -v -n $HOST $PORT << EOF
quote USER $LOGIN
quote PASS $PASSWORD
bin
prompt
get aladin.197901.nc
get analysis.full-arpege.tl149-c24.fa
get glob01.grib
get grid.arome-forecast.guyane0025+0000:00.grib
get grid.arome-forecast.guyane0025+0012:00.grib
get grid.arome-forecast.guyane0025+0024:00.grib
get ic.full-surfex.corsica-02km50.fa
get ICMSHAROM+0022
get pgd.albachir-07km50.SAND.fa
quit
EOF
cd ..
