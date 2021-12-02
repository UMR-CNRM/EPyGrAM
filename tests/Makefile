
NOSE3list     = nosetests-3.7 nosetests3
NOSE3         = $(shell for p in $(NOSE3list) ; do if which $$p >/dev/null 2>&1 ; then echo $$p ; break ; fi ; done)
TEST_ALL      = test_formats.py test_geometries.py test_geometry_methods.py test_spectral.py test_util.py test_combinationsextractions.py
TEST_BASE     = test_formats.py test_geometries.py test_spectral.py test_util.py
NOSE_OPTS     = --verbosity=2 --no-byte-compile
APPTOOLS_DIR  = test_apptools
NB_DIR        = ../epygram/doc_sphinx/source/gallery

.PHONY: all tests tests3 tests_all clean apptools3 apptools

all: tests apptools notebooks_check

# Get data for tests on dedicated repo
get_data:
	. ./get_data.sh

# Run all test with nose
tests: tests3

tests3:
	($(NOSE3) $(TEST_BASE) $(NOSE_OPTS) 2>&1) | tee $@.log

# Complete series
tests_all: tests_all3

tests_all3:
	($(NOSE3) $(TEST_ALL) $(NOSE_OPTS) 2>&1) | tee $@.log

# apptools
apptools: apptools3

apptools3:
	$(MAKE) -C $(APPTOOLS_DIR) all3

# notebooks gallery
notebooks_get_inputs:
	$(MAKE) -C $(NB_DIR) get_inputs

notebooks_check:
	$(MAKE) -C $(NB_DIR) nb_check

clean:
	rm -f *.log
	rm -rf __pycache__
	rm -f *.pyc
