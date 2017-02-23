
/*******************************************************************************/
/* Implémentation parallèle de l'algorithme de bi-orthogonalisation de Lanczos */
/*				 Amira AKLOUL - Maxime KERMARQUER --- M2 CHPS				   */
/*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <mpi.h>
#include <sys/time.h> 
#include <math.h>

double norme_vecteur (double* vecteur, int taille_vecteur){
	double somme = 0;
	double norme;
	int i;

	for (i=0; i<taille_vecteur; i++){
		somme += vecteur[i]*vecteur[i];
	}

	norme = sqrt(somme); 

	return norme;
}

int main(int argc, char const *argv[])
{
	
	return 0;
}