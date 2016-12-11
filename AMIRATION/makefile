all: MATVEC.pgr MATVEC_mpi.pgr

CC= gcc
FLAGS= -g 

MATVEC.pgr: MATVEC.c
	$(CC) $^ $(FLAGS) -fopenmp -o $@

MATVEC_mpi.pgr: MATVEC_mpi.c
	mpicc $^ $(FLAGS) -o $@
exec: MATVEC_mpi.pgr
	mpirun -np 2 ./$^ 4 4 

make clean:
	rm -rf *.pgr
	
