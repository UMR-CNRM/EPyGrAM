EPyGrAM tests
=============

More or less integrated tests for EPyGrAM.

Usage
-----

1. get input data from ftp repo: `make get_data`

2. run basic tests: `make tests`

3. run basic + more advanced tests: `make tests_all`

4. run apptools tests: `make apptools`

5. run notebooks, to check the advanced examples therein work properly: `make notebooks_check`
   (may require to get inputs for notebooks beforehand, with `make notebooks_get_inputs`)

6. clean directory of tests outputs: `make clean`

