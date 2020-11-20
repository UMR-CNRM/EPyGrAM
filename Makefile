# Variables
TEST_DIR	= tests
DOC_DIR		= epygram/doc_sphinx
VGRID_DIR	= site/vgrid/vertical_discretization
SUBDIRS		= tests epygram/doc_sphinx
CLEANDIRS 	= $(SUBDIRS:%=clean-%)

.PHONY: tests doc vgrid mktar pushdev pushrelease clean $(CLEANDIRS)

# TARGETS:
# Run all the test suite using nose
tests:
	$(MAKE) -C $(TEST_DIR) tests
	$(MAKE) -C $(TEST_DIR) apptools

# Build the sphinx documentation
doc:
	$(MAKE) -C $(DOC_DIR)

# Vgrid binary compilation
vgrid:
	$(MAKE) -C $(VGRID_DIR)

# Pushes
pushdev:
	. ./deploy.sh dev

pushrelease:
	. ./deploy.sh

# Tarball for export
tar:
	. ./mktar.sh

# Clean all the directories, then locally
clean: $(CLEANDIRS)
	find . -name "*.pyc"       -print0 | xargs -0r rm
	find . -name "__pycache__" -print0 | xargs -0r rm -r

$(CLEANDIRS):
	$(MAKE) -C $(@:clean-%=%) clean

# Usual target
clobber: clean
