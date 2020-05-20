
PYTHON2 = python2.7
PYTHON3 = python3.7

.PHONY: all all2 all3 cartoplot2 cartoplot3

all: all3

all2: cartoplot2

all3: cartoplot3

cartoplot2:
	. ./test-epy_cartoplot.sh $(PYTHON2)

cartoplot3:
	. ./test-epy_cartoplot.sh $(PYTHON3)
