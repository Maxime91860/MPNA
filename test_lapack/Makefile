

CC = gcc
CCOPT = -O3
LD = gfortran
LDLIBS = -lopenblas

all : col.exe


col.exe : dgeev_col.c
	#$(CC) $^ -o $@ -L/usr/lib/lapack.so -llapack
	$(CC) $^ -o $@ -llapacke




exempledgemm: exempledgemm.o
	$(LD) -o exempledgemm  exempledgemm.o $(LDLIBS)
exempledgemm.o: exempledgemm.c
	$(CC) $(CCOPT) -c exempledgemm.c $(LDLIBS)

clean:
	rm -f exempledgemm  exempledgemm.o


.PHONY: clean