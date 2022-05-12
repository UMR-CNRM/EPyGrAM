# Variables
TEST_DIR	= tests
DOC_DIR		= docs
SUBDIRS		= $(TEST_DIR) $(DOC_DIR)
CLEANDIRS 	= $(SUBDIRS:%=clean-%)

.PHONY: tests doc tar tar_doc tar_tests clear_doc clean $(CLEANDIRS)

# TARGETS:
# Run all the test suite using nose
tests:
	$(MAKE) -C $(TEST_DIR) tests
	$(MAKE) -C $(TEST_DIR) apptools

# Build the sphinx documentation
doc:
	$(MAKE) -C $(DOC_DIR)

# Tarball for export
tar:
	bash ./mktar.sh

# Tarball for export
tar_doc:
	$(MAKE) -C $(DOC_DIR) tar_doc

# Tarball for export
tar_tests:
	$(MAKE) -C $(TEST_DIR) tar_tests

# Clean all the directories, then locally
clean: $(CLEANDIRS)
	find . -name "*.pyc"       -print0 | xargs -0r rm
	find . -name "__pycache__" -print0 | xargs -0r rm -r

$(CLEANDIRS):
	$(MAKE) -C $(@:clean-%=%) clean

clear_doc:
	$(MAKE) -C $(DOC_DIR) clear_doc

# Usual target
clobber: clean
