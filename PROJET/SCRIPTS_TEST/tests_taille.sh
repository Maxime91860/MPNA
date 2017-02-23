#!/bin/bash

for((TAILLE=100;$TAILLE<=2000;TAILLE=TAILLE+100))
	do
	TAILLE_KRYLOV=$(($TAILLE/2))
	echo ""
	echo "./CODE/bi_lanczos_seq.pgr $TAILLE $TAILLE_KRYLOV"
	./CODE/bi_lanczos_seq.pgr $TAILLE $TAILLE_KRYLOV
done

