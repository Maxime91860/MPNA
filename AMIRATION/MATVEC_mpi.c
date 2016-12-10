#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
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
	MPI_Allreduce (&produit_scalaire, &produit_scalaire, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	return resultat;
}




//PRODUIT MATRICE VECTEUR SEQUENTIEL
double *produit_matrice_vecteur(double *matrice, double *vecteur1, int nb_lignes, int nb_colonnes){

	int i;
	double *resultat=(double *)malloc(nb_lignes*sizeof(double));

	for(i=0; i<nb_lignes; i++){
		resultat[i]=produit_scalaire(matrice+(i*nb_colonnes),vecteur1, nb_colonnes);
	}
	return resultat;
}



//MAIN
int main(int argc, char **argv){


	int i, j, nb_lignes, nb_colonnes;

	int size, rank, ideb, ifin, Q, R;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double *vecteur1_complet;
	double *vecteur1;
	double *vecteur2_complet;
	double *vecteur2;
	double *matrice_complete;
	double *matrice;

	struct timeval debut_calcul, fin_calcul, duree_calcul;

	if(argc<3){
		printf("erreur argument\n");
		printf("exemple d'execution: ./executable <nb_lignes> <nb_colonnes>\n"); 
		exit(1);
	}

	nb_lignes=atoi(argv[1]);
	nb_colonnes=atoi(argv[2]);

	Q=nb_lignes/size;
	R=nb_lignes%size;


		if (rank < R){
			ideb = rank * (1+Q);
			ifin = ideb + Q;
		}
		else{
			ideb = R * (1+Q) + ((rank - R) * Q);
			ifin = ideb + Q - 1;
		}


	int nb_lignes_process=ifin-ideb+1;

	if(rank==0){
		vecteur1_complet=(double *)malloc(nb_lignes*sizeof(double));
		vecteur2_complet=(double *)malloc(nb_lignes*sizeof(double));
		matrice_complete=(double *)malloc((nb_lignes*nb_colonnes)*sizeof(double));
	}

	//TOUS LES PROCS ALLOUENT ET INITIALISENT LES VECTEURS ET LA MATRICE

	vecteur1=(double *)malloc(nb_lignes_process*sizeof(double));
	vecteur2=(double *)malloc(nb_lignes_process*sizeof(double));
	matrice=(double *)malloc((nb_lignes_process*nb_colonnes)*sizeof(double));

	for(i=0; i<nb_lignes_process; i++){
		vecteur1[i]=2;
		vecteur2[i]=3;

		for(j=0; j<nb_colonnes; j++){
		matrice[j+(i*nb_colonnes)]=2;		
		}
	}

	MPI_Gather(vecteur1, nb_lignes_process, MPI_DOUBLE, vecteur1_complet, nb_lignes_process, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(vecteur2, nb_lignes_process, MPI_DOUBLE, vecteur2_complet, nb_lignes_process, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if(rank==0){
	printf("vecteur1_complet\n");
	affiche(vecteur1_complet, nb_lignes, nb_lignes);
	printf("vecteur2_complet\n");
	affiche(vecteur2_complet, nb_lignes, nb_lignes);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	

	/*
	printf("vecteur1\n");
	affiche(vecteur1, nb_lignes, nb_lignes);

	printf("vecteur2\n");
	affiche(vecteur2, nb_lignes, nb_lignes);

	printf("matrice\n");
	affiche(matrice, nb_lignes*nb_colonnes, nb_colonnes);
	*/

	//COMPLEXITE EN TEMPS POUR LE PRODUIT SCALAIRE ET LE PRODUIT MATRICE VECTEUR EN SEQUENTIEL



	//COMPLEXITE EN ESPACE MEMOIRE 
	//printf("La complexité en espace memoire : %d\n", (1+ (nb_lignes*3) + (nb_lignes*nb_colonnes)));


	free(vecteur1);
	free(vecteur2);
	free(matrice);
	free(vecteur1_complet);
	free(vecteur2_complet);
	free(matrice_complete);

	MPI_Finalize();
	exit(0);

}
