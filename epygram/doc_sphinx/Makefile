
all: doc clean

.PHONY: all doc clean

doc:
	. ./mk_html_doc.sh

clean:
	find . -name "*.pyc"       -print0 | xargs -0r rm
	find . -name "__pycache__" -print0 | xargs -0r rm -r
	rm -f source/cheatsheet.log
	rm -f source/cheatsheet.aux
	rm -f source/cheatsheet.pdf
