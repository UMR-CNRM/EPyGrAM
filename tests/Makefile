
NOSE3list     = nosetests3 nosetests-3.7
NOSE3         = $(shell for p in $(NOSE3list) ; do if which $$p >/dev/null 2>&1 ; then echo $$p ; break ; fi ; done)
TEST_BASE     = test_formats.py test_geometries.py test_spectral.py test_util.py
TEST_FULL     = $(TEST_BASE) test_geometry_methods.py test_combinationsextractions.py
NOSE_OPTS     = --verbosity=2 --no-byte-compile
APPTOOLS_DIR  = test_apptools
NB_DIR        = ../docs/source/gallery
VERSION       = $(shell head -1 ../VERSION)
TAR_TESTS     = ~/tmp/EPyGrAM-$(VERSION)_tests.tgz

.PHONY: all tests tests_full clean apptools notebooks_check notebooks_get_inputs get_data

all: tests_full apptools notebooks_check

# Get data for tests on dedicated repo
get_data:
	bash ./get_data.sh

# Run tests with nose
tests:
	($(NOSE3) $(TEST_BASE) $(NOSE_OPTS) 2>&1) | tee $@.log

tests_full:
	($(NOSE3) $(TEST_FULL) $(NOSE_OPTS) 2>&1) | tee $@.log

# apptools
apptools:
	$(MAKE) -C $(APPTOOLS_DIR) all

# notebooks gallery
notebooks_get_inputs:
	$(MAKE) -C $(NB_DIR) get_inputs

notebooks_check:
	$(MAKE) -C $(NB_DIR) nb_check

tar_tests: clean
	tar -czf $(TAR_TESTS) *

clean:
	rm -f *.log
	rm -rf __pycache__
	rm -f *.pyc
