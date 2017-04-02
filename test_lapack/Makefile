

CC = gcc
CCOPT = -O3
LD = gfortran
LDLIBS = -lopenblas

exempledgemm: exempledgemm.o
	$(LD) -o exempledgemm  exempledgemm.o $(LDLIBS)
exempledgemm.o: exempledgemm.c
	$(CC) $(CCOPT) -c exempledgemm.c

clean:
	rm -f exempledgemm  exempledgemm.o


.PHONY: clean