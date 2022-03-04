#!/usr/bin/bash

HOST=ftp.umr-cnrm.fr
PORT=21
LOGIN=epygram  # read-only anonymous user
PASSWORD=EpygraM  # with anonymous password

ftp -i -v -n $HOST $PORT << EOF
quote USER $LOGIN
quote PASS $PASSWORD
bin
prompt
cd tests
get data.tgz
quit
EOF

tar -xf data.tgz
rm -f data.tgz

# to put data:
#tar -czf data.tgz data
#ftp $HOST
#cd tests
#put data.tgz
#rm -f data.tgz
