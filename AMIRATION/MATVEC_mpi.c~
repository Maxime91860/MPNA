#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <sys/time.h>



//AFFICHAGE: fonction pour afficher les données
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

//PRODUIT SCALAIRE SEQUENTIEL: fonction du produit scalaire en sequentiel qui retourne un double
double produit_scalaire_sequentiel(double *vecteur1, double *vecteur2, int taille){

	int i;
	double resultat=0;

	for(i=0; i<taille; i++){
		resultat += vecteur1[i]*vecteur2[i];
	}
	
	return resultat;
}

//PRODUIT SCALAIRE PARALLEL: fonction du produit scalaire en parallel qui retourne un double / chaque processus calcul sa partie puis on fait un MPI_Allreduce et on somme tous les elements des vecteurs résultat
double produit_scalaire_parallel(double *vecteur1, double *vecteur2, int taille){

	int i;
	double resultat=0;

	for(i=0; i<taille; i++){
		resultat += vecteur1[i]*vecteur2[i];
	}
	MPI_Allreduce (&resultat, &resultat, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	return resultat;
}



//PRODUIT MATRICE VECTEUR SEQUENTIEL: fonction du produit matrice vecteur en sequentiel qui retourne un double
double *produit_matrice_vecteur_sequentiel(double *matrice, double *vecteur1, int nb_lignes, int nb_colonnes){

	int i;
	double *resultat=(double *)malloc(nb_lignes*sizeof(double));

	for(i=0; i<nb_lignes; i++){
		resultat[i]=produit_scalaire_sequentiel(matrice+(i*nb_colonnes),vecteur1, nb_colonnes);
	}
	return resultat;
	
	}

//PRODUIT MATRICE VECTEUR PARALLEL: chaque processus calcule sa partie du produit matrice vecteur, en sachant que la matrice complete est partitionnée en blocs de lignes dédiés pour chaque processus / on rassemble tous les vecteurs résultat de chaque processus dans le vecteur MATVEC_PARAL avec la fontion MPI_Gather 
void produit_matrice_vecteur_parallel(double *matrice, double *vecteur1, int nb_lignes, int nb_colonnes, int nb_lignes_process, double *MATVEC_PARAL){
					
	int i;
	double *resultat=(double *)malloc(nb_lignes_process*sizeof(double));

	for(i=0; i<nb_lignes_process; i++){
		resultat[i]=produit_scalaire_sequentiel(matrice+(i*nb_colonnes),vecteur1, nb_colonnes);
	}
	MPI_Gather(resultat, nb_lignes_process, MPI_DOUBLE, MATVEC_PARAL, nb_lignes_process, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	free(resultat);
}

//MAIN
int main(int argc, char **argv){

	int size, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//DECLARATION DES VARIABLES
	int i, j, nb_lignes, nb_colonnes;
	int ideb, ifin, Q, R;
	double *vecteur1_complet;
	double *vecteur1;
	double *vecteur2_complet;
	double *vecteur2;
	double *matrice_complete;
	double *matrice;
	double nb_octets=0;
	struct timeval debut_calcul, fin_calcul, duree_calcul;

	//USAGE EXECUTION	
	if(argc<3){
		printf("Erreur argument\n");
		printf("Pour executer le programme en sequentiel: ./executable <nb_lignes> <nb_colonnes>\n"); 
		printf("Pour executer le programme en parallel: mpirun -np <nombre de processus> ./executable <nb_lignes> <nb_colonnes>\n"); 
		exit(1);
	}

	//ARGUMENTS
	nb_lignes=atoi(argv[1]);
	nb_colonnes=atoi(argv[2]);

	//CALCUL DU DEBUT ET DE FIN DE BLOC POUR CHAQUE PROCESSUS
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

	//NOMBRE DE LIGNES DEDIÉES A CHAQUE PROCESSUS
	int nb_lignes_process=ifin-ideb+1;

	//LE RANG 0 ALLOUE DE L'ESPACE MEMOIRE POUR LES VECTEURS COMPLETS ET LA MATRICE COMPLETE
	if(rank==0){
		vecteur1_complet=(double *)malloc(nb_lignes*sizeof(double));
		vecteur2_complet=(double *)malloc(nb_lignes*sizeof(double));
		nb_octets+= 2*nb_lignes*sizeof(double);
		matrice_complete=(double *)malloc((nb_lignes*nb_colonnes)*sizeof(double));
		nb_octets+=(nb_lignes*nb_colonnes)*sizeof(double);
	}

	//TOUS LES PROCS ALLOUENT ET INITIALISENT LEURS PARTIES DE VECTEURS ET LEURS PARTIES DE LA MATRICE
	vecteur1=(double *)malloc(nb_lignes_process*sizeof(double));
	vecteur2=(double *)malloc(nb_lignes_process*sizeof(double));
	nb_octets+= 2*nb_lignes_process*sizeof(double);
	matrice=(double *)malloc((nb_lignes_process*nb_colonnes)*sizeof(double));
	nb_octets+= (nb_lignes_process*nb_colonnes)*sizeof(double);

	for(i=0; i<nb_lignes_process; i++){
		vecteur1[i]= 100;
		vecteur2[i]= 200;

		for(j=0; j<nb_colonnes; j++){
		matrice[j+(i*nb_colonnes)]= 300;		
		}
	}

	//REGROUPEMENT DE TOUTES LES PARTIES QUE DETIENT CHAQUE PROCESSUS DANS DES VECTEURS COMPLETS ET MATRICE COMPLETE DANS LE PROCESSUS DE RANG 0
	MPI_Gather(vecteur1, nb_lignes_process, MPI_DOUBLE, vecteur1_complet, nb_lignes_process, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(vecteur2, nb_lignes_process, MPI_DOUBLE, vecteur2_complet, nb_lignes_process, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(matrice, nb_lignes_process*nb_colonnes, MPI_DOUBLE, matrice_complete, nb_lignes_process*nb_colonnes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	//LE PROCESSUS DE RANG 0 AFFICHE LES DEUX VECTEURS COMPLETS ET LA MATRICE COMPLETE RECCUPÉRÉS DES AUTRES PROCESSUS AVEC LE MPI_Gather
	/*if(rank==0){
	printf("vecteur1_complet\n");
	//affiche(vecteur1_complet, nb_lignes, nb_lignes);
	printf("vecteur2_complet\n");
	//affiche(vecteur2_complet, nb_lignes, nb_lignes);
	printf("la matrice complete\n");
	//affiche(matrice_complete, nb_lignes*nb_colonnes, nb_colonnes);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	*/

//CHAQUE PROCESSUS AFFICHE LES PARTIES DE VECTEURS ET MATRICE QU'IL POSSEDE
	/*for(i=0; i<size; i++){
	 	
	 	printf("Rang : %d debut : %d fin : %d \n",rank,ideb,ifin);
	 	printf("\tvecteur1\n");
	 	affiche(vecteur1, nb_lignes_process, nb_lignes_process);
	 	printf("\tvecteur2\n");
	 	affiche(vecteur2, nb_lignes_process, nb_lignes_process);		
		printf("\tmatrice\n");
		affiche(matrice, nb_lignes_process*nb_colonnes, nb_lignes_process); 
	 }
		MPI_Barrier(MPI_COMM_WORLD);
	*/

	double t1, t2;
	double resultat_parallel, resultat_sequentiel;


	/////////PRODUIT SCALAIRE SEQUENTIEL ET PARALLEL//////////

	//COMPLEXITE EN TEMPS POUR LE PRODUIT SCALAIRE EN PARALLEL
	double temps_produit_matrice_vecteur_parallel;
	double temps_produit_matrice_vecteur_sequentiel;
	double temps_produit_scalaire_sequentiel;
	double temps_produit_scalaire_parallel;
	t1=MPI_Wtime();
	for (i = 0; i < 1000; ++i){
	resultat_parallel=produit_scalaire_parallel(vecteur1, vecteur2, nb_lignes_process);	
	}
	t2=MPI_Wtime();
	temps_produit_scalaire_parallel=t2-t1;

	
	if(rank==0){
		t1=MPI_Wtime();	
		for (i = 0; i < 1000; ++i){		
		resultat_sequentiel=produit_scalaire_sequentiel(vecteur1_complet, vecteur2_complet, nb_lignes);
		}
		t2=MPI_Wtime();
		temps_produit_scalaire_sequentiel=t2-t1;
		printf("\tResultat_sequentiel_produit_scalaire=%f\n",resultat_sequentiel);
		printf("le temps d'execution du produit scalaire en sequentiel est: %f\n", temps_produit_scalaire_sequentiel);
		printf("\tResultat_parallel_produit_scalaire= %f\n",resultat_parallel);
		printf("le temps d'execution du produit scalaire en parallel est: %f\n", temps_produit_scalaire_parallel);
	}
	

	//////////PRODUIT MATRICE VECTEUR SEQUENTIEL ET PARALLEL//////////

	if(rank != 0){
	vecteur1_complet=(double *)malloc(nb_lignes*sizeof(double));
	}
	MPI_Bcast (vecteur1_complet, nb_lignes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	nb_octets+=(size-1)*nb_lignes*sizeof(double);
	MPI_Barrier(MPI_COMM_WORLD);

	double *pmv;
	double *resultat_pmv_sequentiel;
	if(rank ==0){	
	pmv= (double *)malloc(nb_lignes*sizeof(double));	
	resultat_pmv_sequentiel= (double *)malloc(nb_lignes*sizeof(double));
	nb_octets+=nb_lignes*sizeof(double);
	}

	//COMPLEXITE EN TEMPS POUR LE PRODUIT MATRICE VECTEUR EN PARALLEL
	t1=MPI_Wtime();
	for (i = 0; i < 1000; ++i){
	produit_matrice_vecteur_parallel(matrice, vecteur1_complet, nb_lignes, nb_colonnes, nb_lignes_process, pmv);
	}
	t2=MPI_Wtime();
	temps_produit_matrice_vecteur_parallel=t2-t1;

	if(rank==0){
		t1=MPI_Wtime();
		for (i = 0; i < 1000; ++i){
		resultat_pmv_sequentiel=produit_matrice_vecteur_sequentiel(matrice_complete, vecteur1_complet, nb_lignes, nb_colonnes);
		}
		t2=MPI_Wtime();	
		temps_produit_matrice_vecteur_sequentiel=t2-t1;
		//printf("\tResultat_sequentiel_produit_matrice_vecteur1_complet = \n");
		//affiche(resultat_pmv_sequentiel, nb_lignes, nb_lignes);
		printf("le temps d'execution du produit matrice vecteur en sequentiel est: %f\n", temps_produit_matrice_vecteur_sequentiel);
		//printf("\tResultat_parallel_produit_matrice_vecteur1_complet = \n");	
		//affiche(pmv, nb_lignes, nb_lignes);
		printf("le temps d'execution du produit matrice vecteur en parallel est: %f\n", temps_produit_matrice_vecteur_parallel);
	}


	//COMPLEXITE EN ESPACE MEMOIRE 
	/*double nb_total_octets = 0;

	MPI_Reduce (&nb_octets, &nb_total_octets, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	*/

	if(rank ==0){
		//printf("La complexité en espace memoire : %d\n", (1+ (nb_lignes*3) + (nb_lignes*nb_colonnes)));
		printf("Le nombre d'octets total alloué : %f octets.\n", nb_octets);	
		
	FILE* out = fopen("resultat.txt","a+");
		//NB_PROCESSUS - TAILLE MATRICE -TEMPS SEQUENTIEL PRODUIT SCALAIRE -TEMPS PARALLEL PRODUIT SCALAIRE - TEMPS SEQUENTIEL PRODUIT MATRICE VECTEUR - TEMPS PARALLEL PRODUIT MATRICE VECTEUR
		fprintf(out, "%d %d %f %f %f %f\n", size, nb_lignes, temps_produit_scalaire_sequentiel, temps_produit_scalaire_parallel, temps_produit_matrice_vecteur_sequentiel, temps_produit_matrice_vecteur_parallel);
	
	}

	free(vecteur1);
	free(vecteur2);
	free(matrice);
	
	if(rank ==0){	
		free(vecteur1_complet);
		free(vecteur2_complet);
		free(matrice_complete);
		free(pmv);	
		free(resultat_pmv_sequentiel);
	}

	MPI_Finalize();
	exit(0);

}
