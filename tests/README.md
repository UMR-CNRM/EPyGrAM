EPyGrAM tests
=============

More or less integrated tests for EPyGrAM.

Usage
-----

1. the first time, get input data from ftp repo: `make get_data`

2. run tests:
   * basic tests: `make tests`
   * basic + more advanced tests: `make tests_full`
   * apptools tests: `make apptools`
   * run notebooks, to check the advanced examples therein work properly: `make notebooks_check`
     (may require to get inputs for notebooks beforehand, the first time, with `make notebooks_get_inputs`)
   * or them all at once: `make all`

3. clean directory of tests outputs: `make clean`

