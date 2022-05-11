# Variables
TEST_DIR	= tests
DOC_DIR		= epygram/doc_sphinx
SUBDIRS		= tests epygram/doc_sphinx
CLEANDIRS 	= $(SUBDIRS:%=clean-%)

.PHONY: notebooks4doc notebooks_clean tests doc mktar clean $(CLEANDIRS)

# TARGETS:
# Run all the test suite using nose
tests:
	$(MAKE) -C $(TEST_DIR) tests
	$(MAKE) -C $(TEST_DIR) apptools

# Build the sphinx documentation
doc:
	$(MAKE) -C $(DOC_DIR)

notebooks4doc:
	$(MAKE) -C $(DOC_DIR) nb4doc

notebooks_clean:
	$(MAKE) -C $(DOC_DIR) clean

# Tarball for export
tar:
	bash ./mktar.sh

# Tarball for export
tar_doc:
	bash ./mktar_doc.sh

# Tarball for export
tar_tests:
	bash ./mktar_tests.sh

# Clean all the directories, then locally
clean: $(CLEANDIRS)
	find . -name "*.pyc"       -print0 | xargs -0r rm
	find . -name "__pycache__" -print0 | xargs -0r rm -r

$(CLEANDIRS):
	$(MAKE) -C $(@:clean-%=%) clean

# Usual target
clobber: clean
