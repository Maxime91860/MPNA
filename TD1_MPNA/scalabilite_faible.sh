#!/bin/bash
NBTACHE=4
rm -rf resultat.txt

for((TAILLE=100;$TAILLE<1000;TAILLE=TAILLE+100))
	do
	echo "mpirun -np $NBTACHE ./td1_mpna_mpi.pgr $TAILLE"
	mpirun -np $NBTACHE  ./td1_mpna_mpi.pgr $TAILLE
done

for((TAILLE=1000;$TAILLE<=15000;TAILLE=TAILLE+1000))
	do
	echo "mpirun -np $NBTACHE ./td1_mpna_mpi.pgr $TAILLE"
	mpirun -np $NBTACHE  ./td1_mpna_mpi.pgr $TAILLE
done