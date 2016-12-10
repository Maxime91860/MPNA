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
double produit_scalaire_sequentiel(double *vecteur1, double *vecteur2, int taille){

	int i;
	double resultat=0;

	for(i=0; i<taille; i++){
		resultat += vecteur1[i]*vecteur2[i];
	}
	
	return resultat;
}

double produit_scalaire_parallel(double *vecteur1, double *vecteur2, int taille){

	int i;
	double resultat=0;

	for(i=0; i<taille; i++){
		resultat += vecteur1[i]*vecteur2[i];
	}
	MPI_Allreduce (&resultat, &resultat, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	return resultat;
}



//PRODUIT MATRICE VECTEUR SEQUENTIEL
double *produit_matrice_vecteur(double *matrice, double *vecteur1, int nb_lignes, int nb_colonnes){

	int i;
	double *resultat=(double *)malloc(nb_lignes*sizeof(double));

	for(i=0; i<nb_lignes; i++){
		resultat[i]=produit_scalaire_parallel(matrice+(i*nb_colonnes),vecteur1, nb_colonnes);
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
	MPI_Gather(matrice, nb_lignes_process*nb_colonnes, MPI_DOUBLE, matrice_complete, nb_lignes_process*nb_colonnes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	if(rank==0){
	printf("vecteur1_complet\n");
	affiche(vecteur1_complet, nb_lignes, nb_lignes);
	printf("vecteur2_complet\n");
	affiche(vecteur2_complet, nb_lignes, nb_lignes);
	printf("la matrice complete\n");
	affiche(matrice_complete, nb_lignes*nb_colonnes, nb_colonnes);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	int t1, t2;
	double resultat_parallel, resultat_sequentiel;

	t1=MPI_Wtime();
	resultat_parallel=produit_scalaire_parallel(vecteur1, vecteur2, nb_lignes_process);	
	t2=MPI_Wtime();

	
	if(rank==0){
		resultat_sequentiel=produit_scalaire_sequentiel(vecteur1_complet, vecteur2_complet, nb_lignes);
		printf("\tResultat sequentiel=%f\n",resultat_sequentiel);
		printf("\tResultat_parallel= %f\n",resultat_parallel);
		printf("le temps d'execution du produit scalaire en parallel est: %d\n", t2-t1);
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



	//COMPLEXITE EN ESPACE MEMOIRE 
	//printf("La complexité en espace memoire : %d\n", (1+ (nb_lignes*3) + (nb_lignes*nb_colonnes)));
	
	

	free(vecteur1);
	free(vecteur2);
	free(matrice);
	
	if(rank ==0){	
		free(vecteur1_complet);
		free(vecteur2_complet);
		free(matrice_complete);
	}

	MPI_Finalize();
	exit(0);

}
