FC=gfortran
FLAGS="-fPIC"
SHARED_FLAGS="-shared"
EXEC=mkvgrid.x
SO=atmosta.so

all: $(EXEC) $(SO)

mkvgrid.o: mkvgrid.F90
	$(FC) $(FLAGS) -o mkvgrid.o -c mkvgrid.F90

altsta.o: altsta.F90
	$(FC) $(FLAGS) -o altsta.o -c altsta.F90

presta.o: presta.F90 altsta.o
	$(FC) $(FLAGS) -o presta.o -c presta.F90

posname_nodep.o: posname_nodep.F90
	$(FC) $(FLAGS) -o posname_nodep.o -c posname_nodep.F90

mkvgrid.x: mkvgrid.o posname_nodep.o altsta.o presta.o
	$(FC) $(FLAGS) -o mkvgrid.x posname_nodep.o altsta.o presta.o mkvgrid.o

atmosta.so: altsta.o presta.o
	$(FC) $(FLAGS) $(SHARED_FLAGS) -o atmosta.so altsta.o presta.o

clean:
	rm -rf *.o

mrproper: clean cleanout
	rm -rf $(EXEC) $(SO)

cleanout:
	rm -rf *.namvv1 *.infosup *.mapping *.latex *.html
