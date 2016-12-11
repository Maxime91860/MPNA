all: MATVEC.pgr MATVEC_mpi.pgr

CC= gcc
FLAGS= -g 

MATVEC.pgr:MATVEC.c
	$(CC) $^ $(FLAGS) -fopenmp -o $@

MATVEC_mpi.pgr: MATVEC_mpi.c
	mpicc $^ $(FLAGS) -o $@
exec:MATVEC_mpi.pgr
	mpirun -np 4 ./$^ 12 12 

make clean:
	rm -rf *.pgr
	
