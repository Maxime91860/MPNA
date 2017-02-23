#!/bin/bash

for((TAILLE=10;$TAILLE<=100;TAILLE=TAILLE+10))
	do
	TAILLE_KRYLOV=$TAILLE/2
	echo ""
	echo "../CODE/bi_lanczos_seq.pgr $TAILLE $TAILLE_KRYLOV"
	../CODE/bi_lanczos_seq.pgr $TAILLE $TAILLE_KRYLOV
done

