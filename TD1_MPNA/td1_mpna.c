
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <omp.h>
#include <sys/time.h> 


void affiche(double* tab, int N, int width){
	int i;
	printf("(");
	for(i=0; i<N; i++){
		if((i+1)%width == 0){
			printf("%f)\n",tab[i]);
			if(i != (N-1)){
				printf("(");
			}
		}
		else{
			printf("%f ",tab[i]);
		}
	}
}

double produit_scalaire (double* vecteur1, double* vecteur2, int taille_vecteur){
	double produit_scalaire;
	int i;
	produit_scalaire = 0;

	for(i=0; i<taille_vecteur; i++){
		produit_scalaire += vecteur1[i] * vecteur2[i];
	}

	return produit_scalaire;
}


double produit_scalaire_parralle (double* vecteur1, double* vecteur2, int taille_vecteur){
	double produit_scalaire;
	int i;
	produit_scalaire = 0;


	for(i=0; i<taille_vecteur; i++){
		produit_scalaire += vecteur1[i] * vecteur2[i];
	}

	return produit_scalaire;
}


double* produit_matrice_vecteur (double* matrice, double* vecteur, int nb_ligne, int nb_col){
	double* resultat;
	int i;

	resultat = (double*) malloc(nb_ligne*sizeof(double));

	for(i=0; i<nb_ligne; i++){
		resultat[i] = produit_scalaire(matrice+(i*nb_col), vecteur, nb_col);
	}

	return resultat;
}

int main (int argc, char** argv){

	if(argc < 3){
		fprintf(stderr, "Erreur argument : Usage %s <int : nbligne> <int : nbcol>\n",argv[0]);
		exit(-1);
	}

	int nb_ligne, nb_col;
	double* vecteur1;
	double* vecteur2;
	double*  matrice;
	double*  matrice_vecteur;
	int i,j;

	struct timeval debut_calcul, fin_calcul, duree_calcul;

	nb_ligne = atoi(argv[1]);
	nb_col = atoi(argv[2]);

	vecteur1 = (double*) malloc(nb_ligne*sizeof(double));
	vecteur2 = (double*) malloc(nb_ligne*sizeof(double));
	matrice  = (double*) malloc(nb_ligne*nb_col*sizeof(double));

	//Initialisation
	for(i=0; i<nb_ligne; i++){
		vecteur1[i] = (double) (rand()%100);
		vecteur2[i] = (double) (rand()%100);
		for(j=0; j<nb_col; j++){
			matrice[j + i*nb_col] = (double) i;
		}
	}

	printf("Vecteur 1\n");
	// affiche(vecteur1, nb_ligne, nb_ligne);
	printf("Vecteur 2\n");
	// affiche(vecteur2, nb_ligne, nb_ligne);

	printf("Matrice\n");
	// affiche(matrice, nb_ligne*nb_col, nb_col);

	// printf("\nProduit scalaire de vecteur1 par vecteur2 : %f\n", produit_scalaire(vecteur1, vecteur2, nb_ligne));
	gettimeofday(&debut_calcul, NULL);
	produit_scalaire(vecteur1, vecteur2, nb_ligne);
	gettimeofday(&fin_calcul, NULL);
	timersub(&fin_calcul, &debut_calcul, &duree_calcul);
	printf("Temps produit scalaire sequentiel : %ld.%ld s\n",(long int)duree_calcul.tv_sec,(long int)duree_calcul.tv_usec);


	// printf("\nProduit de matrice par vecteur1 : \n");
	// matrice_vecteur = produit_matrice_vecteur(matrice, vecteur1, nb_ligne, nb_col);
	// affiche(matrice_vecteur, nb_ligne, 1);
	gettimeofday(&debut_calcul, NULL);
	produit_matrice_vecteur(matrice, vecteur1, nb_ligne, nb_col);
	gettimeofday(&fin_calcul, NULL);
	timersub(&fin_calcul, &debut_calcul, &duree_calcul);
	printf("Temps matrice vecteur sequentiel : %ld.%ld s\n",(long int)duree_calcul.tv_sec,(long int)duree_calcul.tv_usec);


	free(vecteur1);
	free(vecteur2);
	free(matrice);

	exit(0);
}