
all : td1_mpna.pgr

td1_mpna.pgr : td1_mpna.c
	gcc $^ -g -fopenmp -Wall -o $@

clean : 
	rm -rf *.pgr