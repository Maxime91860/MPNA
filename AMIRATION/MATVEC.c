#include <stdio.h>
#include <stdlib.h>
//#include <mpi.h>
#include <omp.h>
#include <sys/time.h>



//AFFICHAGE
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

//PRODUIT SCALAIRE SEQUENTIEL
double produit_scalaire(double *vecteur1, double *vecteur2, int taille){

int i;
double resultat=0;

for(i=0; i<taille; i++){
	resultat += vecteur1[i]*vecteur2[i];
}
return resultat;
}

//PRODUIT SCALAIRE PARALLEL
double produit_scalaire_parallel(double *vecteur1, double *vecteur2, int taille){

int i;
double resultat=0;

#pragma omp parallel private(i) shared(vecteur1, vecteur2, taille)
#pragma omp for
for(i=0; i<taille; i++){
	resultat+=vecteur1[i] * vecteur2[i];
}
return resultat;
}


//PRODUIT MATRICE VECTEUR SEQUENTIEL
double *produit_matrice_vecteur(double *matrice, double *vecteur1, int nb_lignes, int nb_colonnes){

int i;
double *resultat=(double *)malloc(nb_lignes*sizeof(double));

for(i=0; i<nb_lignes; i++){
	resultat[i]+=produit_scalaire(matrice+(i*nb_colonnes),vecteur1, nb_colonnes);
}
return resultat;
}

//PRODUIT MATRICE VECTEUR PARALLEL

double *produit_matrice_vecteur_parallel(double *matrice, double *vecteur1, int nb_lignes, int nb_colonnes){

int i;
double *resultat=(double *)malloc(nb_lignes*nb_colonnes*sizeof(double));

#pragma omp parallel private(i) shared(matrice, vecteur1, nb_lignes, nb_colonnes)
#pragma omp for
for(i=0; i<nb_lignes; i++){
	resultat[i]+=produit_scalaire_parallel(matrice+i*nb_colonnes, vecteur1, nb_colonnes);
}
return resultat;
}



//MAIN
int main(int argc, char **argv){

int i, j, nb_lignes, nb_colonnes;

struct timeval debut_calcul, fin_calcul, duree_calcul;

if(argc<3){
	printf("erreur argument\n");
	printf("exemple d'execution: ./executable <nb_lignes> <nb_colonnes>\n"); 
	exit(1);
}

nb_lignes=atoi(argv[1]);
nb_colonnes=atoi(argv[2]);

double *vecteur1=(double *)malloc(nb_lignes*sizeof(double));
double *vecteur2=(double *)malloc(nb_lignes*sizeof(double));
double *matrice=(double *)malloc((nb_lignes*nb_colonnes)*sizeof(double));

//INITIALISATION DES VECTEURS ET DES MATRICES

printf("Initialisation des données\n");

for(i=0; i<nb_lignes; i++){
	vecteur1[i]=2;
	vecteur2[i]=2;

	for(j=0; j<nb_colonnes; j++){
	matrice[j+(i*nb_colonnes)]=2;		
	}
}

/*
printf("vecteur1\n");
affiche(vecteur1, nb_lignes, nb_lignes);

printf("vecteur2\n");
affiche(vecteur2, nb_lignes, nb_lignes);

printf("matrice\n");
affiche(matrice, nb_lignes*nb_colonnes, nb_colonnes);
*/

//COMPLEXITE EN TEMPS POUR LE PRODUIT SCALAIRE ET LE PRODUIT MATRICE VECTEUR EN SEQUENTIEL

gettimeofday(&debut_calcul, NULL);
produit_scalaire(vecteur1, vecteur2, nb_lignes);
gettimeofday(&fin_calcul, NULL);
timersub(&fin_calcul, &debut_calcul, &duree_calcul);
printf("temps d'execution du produit scalaire en séquentiel: %fs\n", (double) (duree_calcul.tv_sec)+(duree_calcul.tv_usec) /1000000.0);



gettimeofday(&debut_calcul, NULL);
produit_matrice_vecteur(matrice, vecteur1, nb_lignes, nb_colonnes);
gettimeofday(&fin_calcul, NULL);
timersub(&fin_calcul, &debut_calcul, &duree_calcul);
printf("le temps d'execution du produit matrice vecteur en séquentiel: %fs\n", (double) (duree_calcul.tv_sec)+(duree_calcul.tv_usec) / 1000000.0);


//COMPLEXITE EN TEMPS POUR LE PRODUIT SCALAIRE ET LE PRODUIT MATRICE VECTEUR EN PARALLEL

gettimeofday(&debut_calcul, NULL);
produit_scalaire_parallel(vecteur1, vecteur2, nb_lignes);
gettimeofday(&fin_calcul, NULL);
timersub(&fin_calcul, &debut_calcul, &duree_calcul);
printf("temps d'execution du produit scalaire en parallel: %fs\n", (double) (duree_calcul.tv_sec)+(duree_calcul.tv_usec) /1000000.0);



gettimeofday(&debut_calcul, NULL);
produit_matrice_vecteur_parallel(matrice, vecteur1, nb_lignes, nb_colonnes);
gettimeofday(&fin_calcul, NULL);
timersub(&fin_calcul, &debut_calcul, &duree_calcul);
printf("le temps d'execution du produit matrice vecteur en parallel: %fs\n", (double) (duree_calcul.tv_sec)+(duree_calcul.tv_usec) / 1000000.0);


//COMPLEXITE EN ESPACE MEMOIRE 
printf("La complexité en espace memoire : %d\n", (1+ (nb_lignes*3) + (nb_lignes*nb_colonnes)));


free(vecteur1);
free(vecteur2);
free(matrice);

exit(0);

}
