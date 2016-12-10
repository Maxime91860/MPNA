
all : td1_mpna.pgr td1_mpna_mpi.pgr

FLAGS=-g


td1_mpna.pgr : td1_mpna.c
	gcc $^ $(FLAGS) -fopenmp -o $@

td1_mpna_mpi.pgr : td1_mpna_mpi.c
	mpicc $^ $(FLAGS) -o $@

exec : td1_mpna_mpi.pgr 
	mpirun -np 2 ./td1_mpna_mpi.pgr 4 4

clean : 
	rm -rf *.pgr
