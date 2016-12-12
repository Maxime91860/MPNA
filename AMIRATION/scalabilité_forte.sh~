#!/bin/bash
TAILLE=15000
rm -rf resultat.txt

for((NBTACHE=1;$NBTACHE<4;NBTACHE=NBTACHE+1))
	do
	echo "mpirun -np $NBTACHE ./MATVEC_mpi.pgr $TAILLE"
	mpirun -np $NBTACHE  ./MATVEC_mpi.pgr $TAILLE $TAILLE
done

