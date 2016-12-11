#!/bin/bash
TAILLE=15000
rm -rf resultat.txt

for((NBTACHE=1;$NBTACHE<10;NBTACHE=NBTACHE+1))
	do
	echo "mpirun -np $NBTACHE ./td1_mpna_mpi.pgr $TAILLE"
	mpirun -np $NBTACHE  ./td1_mpna_mpi.pgr $TAILLE
done
