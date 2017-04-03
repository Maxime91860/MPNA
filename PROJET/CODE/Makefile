

MPI_PATH=`dirname \` which mpicc \` `
MPI_FLAGS=-I$(MPI_PATH)/../include/mpi -L$(MPI_PATH)/../lib -lmpi
LAPACK_FLAGS=-llapacke
#-lblas 
LIB_FLAGS=-g -Wall -lm
CC=gcc

all : bi_lanczos_seq.pgr #bi_lanczos_par.pgr

bi_lanczos_seq.pgr : bi_lanczos_seq.c
	$(CC) $^ $(LIB_FLAGS) $(LAPACK_FLAGS) -o $@

bi_lanczos_par.pgr : bi_lanczos_par.c
	$(CC) $^ $(LIB_FLAGS) $(MPI_FLAGS) -o $@

clean :
	rm -rf *.pgr