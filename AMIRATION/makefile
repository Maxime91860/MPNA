all: MATVEC.pgr

CC= gcc
FLAGS= -g -Wall

MATVEC.pgr:MATVEC.c
	$(CC) $^ $(FLAGS) -fopenmp -o $@

make clean:
	rm -rf *.pgr
	
