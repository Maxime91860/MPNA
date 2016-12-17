
all : td1_mpna.pgr td1_mpna_mpi.pgr gram_schmidt.pgr	

FLAGS=-g -Wall


td1_mpna.pgr : td1_mpna.c
	gcc $^ $(FLAGS) -fopenmp -o $@

td1_mpna_mpi.pgr : td1_mpna_mpi.c
	mpicc $^ $(FLAGS) -o $@

gram_schmidt.pgr : gram_schmidt.c
	gcc $^ $(FLAGS) -lm -o $@

exec : td1_mpna_mpi.pgr 
	mpirun -np 3 ./td1_mpna_mpi.pgr 3

clean : 
	rm -rf *.pgr
