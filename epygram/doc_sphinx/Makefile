
NB_DIR		= source/gallery

all: doc clean

.PHONY: all doc clean notebooks4doc

doc:
	. ./mk_html_doc.sh

notebooks4doc:
	$(MAKE) -C $(NB_DIR) nb4doc

get_gallery_inputs:
	$(MAKE) -C $(NB_DIR) get_inputs

clean:
	find . -name "*.pyc"       -print0 | xargs -0r rm
	find . -name "__pycache__" -print0 | xargs -0r rm -r
	rm -f source/cheatsheet.log
	rm -f source/cheatsheet.aux
	rm -f source/cheatsheet.pdf

clean_doc:
	rm -rf html

