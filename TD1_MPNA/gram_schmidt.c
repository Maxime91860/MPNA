
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/time.h> 
#include <math.h>


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



double* produit_matrice_vecteur (double* matrice, double* vecteur, int nb_ligne, int nb_col){
	double* resultat;
	int i;

	resultat = (double*) malloc(nb_ligne*sizeof(double));

	for(i=0; i<nb_ligne; i++){
		resultat[i] = produit_scalaire(matrice+(i*nb_col), vecteur, nb_col);
	}

	return resultat;
}

double norme_vecteur (double* vecteur, int taille_vecteur){
	double somme = 0;
	int i;

	for (i=0; i<taille_vecteur; i++){
		somme += vecteur[i]*vecteur[i];
	}

	somme = sqrt(somme); 

	return somme;
}


double* gram_schmidt_modifie (double* A, int nb_vecteurs, int taille_vecteur){

	double* Q = (double *) malloc (nb_vecteurs*taille_vecteur*sizeof(double));
	double* w;
	double* R = (double *) malloc (nb_vecteurs*taille_vecteur*sizeof(double));
	int i, j, k;

	for(i=0; i<nb_vecteurs*taille_vecteur; i++){
		R[i] = 0;
	}

	for (k=0; k<nb_vecteurs; k++){
		w = A+(k*taille_vecteur);
		for(j=0; j<=k-1; j++){
			R[j + k*taille_vecteur] = produit_scalaire(w, Q+j*taille_vecteur, taille_vecteur);
		
			for(i=0; i<taille_vecteur; i++){
				w[i] = w[i] - R[j + k*taille_vecteur] * Q[i + j*taille_vecteur];
			}
		}
		R[k + k*taille_vecteur] = norme_vecteur(w,taille_vecteur);
		for(i=0; i<taille_vecteur; i++){
			Q[i + k*taille_vecteur] = w[i] / R[k + k*taille_vecteur];
		}

		// printf("\n");
		// affiche(R, nb_vecteurs*taille_vecteur, taille_vecteur);
		// printf("\n");
	}
	return Q;
}

double* gram_schmidt (double* A, int nb_vecteurs, int taille_vecteur){

	double* Q = (double *) malloc (nb_vecteurs*taille_vecteur*sizeof(double));
	double* w;
	double* R = (double *) malloc (nb_vecteurs*taille_vecteur*sizeof(double));
	int i, j, k;

	for(i=0; i<nb_vecteurs*taille_vecteur; i++){
		R[i] = 0;
	}

	for (k=0; k<nb_vecteurs; k++){
		w = A+(k*taille_vecteur);
		for(j=0; j<=k-1; j++){
			R[j + k*taille_vecteur] = produit_scalaire(w, Q+j*taille_vecteur, taille_vecteur);
		}
		for(j=0; j<=k-1; j++){
			for(i=0; i<taille_vecteur; i++){
				w[i] = w[i] - R[j + k*taille_vecteur] * Q[i + j*taille_vecteur];
			}
		}
		R[k + k*taille_vecteur] = norme_vecteur(w,taille_vecteur);
		for(i=0; i<taille_vecteur; i++){
			Q[i + k*taille_vecteur] = w[i] / R[k + k*taille_vecteur];
		}

		// printf("\n");
		// affiche(R, nb_vecteurs*taille_vecteur, taille_vecteur);
		// printf("\n");
	}
	return Q;
}



int main (int argc, char** argv){

	if(argc < 2){
		fprintf(stderr, "Erreur argument : Usage %s <int : Taille>\n",argv[0]);
		exit(-1);
	}

	int nb_vecteurs;
	int taille_vecteur;
	double* matrice_A;
	double* resultat;
	double* resultat_modifie;
	int i, j;

	nb_vecteurs = atoi(argv[1]);
	taille_vecteur = atoi(argv[1]);

	//Allocation 
	matrice_A = (double *) malloc (nb_vecteurs*taille_vecteur*sizeof(double));

	//Initialisation
	for(i=0; i<nb_vecteurs*taille_vecteur; i++){
		matrice_A[i] = random()/(double)RAND_MAX;
	}

	// affiche(matrice_A, nb_vecteurs*taille_vecteur, taille_vecteur);

	resultat = gram_schmidt(matrice_A, nb_vecteurs, taille_vecteur);
	resultat_modifie = gram_schmidt_modifie(matrice_A, nb_vecteurs, taille_vecteur);

	// printf("\n");

	// affiche(resultat, nb_vecteurs*taille_vecteur, taille_vecteur);
	// printf("\n");

	double erreur_precision = 0;
	for(i=0; i<nb_vecteurs; i++){
		for(j=0+i; j<nb_vecteurs; j++){
			if(i != j){
				erreur_precision += fabs(produit_scalaire(resultat + i*taille_vecteur, 
														  resultat + j*taille_vecteur, taille_vecteur));
			}
		}
	}

	double erreur_precision_modfie = 0;
	for(i=0; i<nb_vecteurs; i++){
		for(j=0+i; j<nb_vecteurs; j++){
			if(i != j){
				erreur_precision_modfie += fabs(produit_scalaire(resultat_modifie + i*taille_vecteur,
																 resultat_modifie + j*taille_vecteur, taille_vecteur));
			}
		}
	}



	printf("Erreur de precision Gram-Schmidt: %e\n",erreur_precision);
	printf("Erreur de precision Gram-Schmidt modifié: %e\n",erreur_precision_modfie);

	FILE* fichier_sortie = fopen("gram_schmidt.txt","a+");
	fprintf(fichier_sortie, "%d %e %e\n", nb_vecteurs, erreur_precision, erreur_precision_modfie);


	free(resultat);
	free(matrice_A);


	exit(0);
}


