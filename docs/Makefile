
NB_DIR		= source/gallery

all: doc clean

.PHONY: all doc clean get_gallery_inputs notebooks_clean clear_doc

doc:
	bash ./mk_html_doc.sh

get_gallery_inputs:
	$(MAKE) -C $(NB_DIR) get_inputs

notebooks_clean:
	$(MAKE) -C $(NB_DIR) clean

clean: notebooks_clean
	rm -f source/cheatsheet.log
	rm -f source/cheatsheet.aux
	rm -f source/cheatsheet.pdf

clear_doc:
	rm -rf html
	rm -rf build

