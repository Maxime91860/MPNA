
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <mpi.h>
#include <sys/time.h> 


// struct decoupage_mpi {
// 	int debut;
// 	int fin;
// };
// struct decoupage_mpi decoupage = (struct decoupage_mpi) malloc (sizeof(struct decoupage_mpi));

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

void produit_matrice_vecteur (double* matrice, double* vecteur, int nb_ligne, int nb_col, double* PMV){
	int i;

	for(i=0; i<nb_ligne; i++){
		PMV[i] = produit_scalaire(matrice+(i*nb_col), vecteur, nb_col);
	}

}

void produit_scalaire_parallele (double* vecteur1, double* vecteur2, int taille_vecteur, double* produit_scalaire){
	int i;
	*produit_scalaire = 0;

	//Chaque processus va faire une partie des calculs
	for(i=0; i<taille_vecteur; i++){
		*produit_scalaire += vecteur1[i] *  vecteur2[i];
	}

	MPI_Allreduce (produit_scalaire, produit_scalaire, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void produit_matrice_vecteur_parallele (double* matrice, double* vecteur, int nb_ligne, int nbligne_process, int nb_col, double* PMV){
	int i;
	double* resultat;
	int id_process;
	MPI_Comm_rank(MPI_COMM_WORLD,&id_process);

	resultat = (double*) malloc(nbligne_process*sizeof(double));

	for(i=0; i<nbligne_process; i++){
		resultat[i] = produit_scalaire(matrice+(i*nb_col), vecteur, nb_col);
	}


	MPI_Gather(resultat, nbligne_process, MPI_DOUBLE, PMV, nbligne_process, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	free(resultat);
}

int main (int argc, char** argv){
     
	if(argc < 3){
		fprintf(stderr, "Erreur argument : Usage %s <int : nbligne> <int : nbcol>\n",argv[0]);
		exit(-1);
	}
 
	int nb_ligne, nb_col;
	double* vecteur1;
	double* vecteur1_complet;
	double* vecteur2;
	double* vecteur2_complet;
	double* matrice;
	double* matrice_complete;
	int i,j;

	int size, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	nb_ligne = atoi(argv[1]);
	nb_col = atoi(argv[2]);


	//Le processus 0 va faire les calculs en sequentiel pour pouvoir comparer ensuite les
	//resultats sequentiels et paralleles.
	if(rank==0){
		vecteur1_complet = (double*) malloc(nb_ligne*sizeof(double));
		vecteur2_complet = (double*) malloc(nb_ligne*sizeof(double));
		matrice_complete = (double*) malloc(nb_ligne*nb_col*sizeof(double));
	}


	//Calcul des decoupages des données 
	int debut, fin;
	int quotient, reste;
	quotient = nb_ligne / size;
	reste = nb_ligne % size;

	if (rank < reste){
		debut = rank * (1+quotient);
		fin = debut + quotient;
	}
	else{
		debut = reste * (1+quotient) + ((rank - reste) * quotient);
		fin = debut + quotient - 1;
	}

	int nbligne_process = fin-debut+1;


	//Allocation par tous les processus
	vecteur1 = (double*) malloc(nbligne_process*sizeof(double));
	vecteur2 = (double*) malloc(nbligne_process*sizeof(double));
	matrice  = (double*) malloc(nbligne_process*nb_col*sizeof(double));

	//Initialisation par tous les processus
	for(i=0; i<nbligne_process; i++){
		vecteur1[i] = 2; //(double) (rand()%100);
		vecteur2[i] = 3; //(double) (rand()%100);
		for(j=0; j<nb_col; j++){
			if((i+debut)==j){
				matrice[j + i*nb_col] = rank+1;	
			}
			else{
				matrice[j + i*nb_col] = 0;		
			}
			
		}
	}


	//Affichage des des données de chaque processus
	// for(i=0; i<size; i++){
	// 	if(i==rank){
	// 		printf("Rang : %d\n",rank);
	// 		printf("\tdebut = %d, fin = %d\n",debut,fin);
	// 		printf("\tVecteur1\n");
	// 		affiche(vecteur1, nbligne_process , nbligne_process);
	// 		printf("\tVecteur2\n");
	// 		affiche(vecteur2, nbligne_process, nbligne_process);
	// 		affiche(matrice, nbligne_process*nb_col, nb_col);
	// 		printf("\n");
	// 	}
	// 	MPI_Barrier(MPI_COMM_WORLD);
	// }
	// MPI_Barrier(MPI_COMM_WORLD);


	//Le processus 0 récupère les données 
	MPI_Gather(vecteur1, nbligne_process, MPI_DOUBLE, vecteur1_complet, nbligne_process, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(vecteur2, nbligne_process, MPI_DOUBLE, vecteur2_complet, nbligne_process, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(matrice, nbligne_process*nb_col, MPI_DOUBLE, matrice_complete, nbligne_process*nb_col, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	//Affichage des données récupérées par le processus 0
	if(rank == 0){
		printf("Vecteur1 complet\n");
		//affiche(vecteur1_complet, nb_ligne, nb_ligne);
		printf("Vecteur2 complet\n");
		//affiche(vecteur2_complet, nb_ligne, nb_ligne);
		printf("Matice complete\n");
		//affiche(matrice_complete, nb_ligne*nb_col, nb_col);
	}
	MPI_Barrier(MPI_COMM_WORLD);


  	//---------------------------------------------------------------------------------------------//

	double t1, t2;

	/* ******************************* */
	//	 	 PRODUIT SCALAIRE 			/
	/* ******************************* */
	double resultat;
	double resultat_seq = 0;

	t1 = MPI_Wtime();
	for (i = 0; i < 1000; ++i)
	{
		produit_scalaire_parallele(vecteur1, vecteur2, nbligne_process, &resultat);
	}
	t2 = MPI_Wtime();

	if(rank == 0){
		printf("\n --- Produit scalaire vecteur1 * vecteur2 --- \n");
		printf("\nTemps produit scalaire parallèle : %f s\n", t2-t1);
		t1 = MPI_Wtime();
		for (i = 0; i < 1000; ++i)
		{
			produit_scalaire(vecteur1_complet , vecteur2_complet, nb_ligne);
		}
		t2 = MPI_Wtime();
		resultat_seq = produit_scalaire(vecteur1_complet , vecteur2_complet, nb_ligne);
		printf("\tResultat parallèle = %f\n",resultat);
		printf("Temps produit scalaire séquentiel : %f s\n", t2-t1);
		printf("\tResultat séquentiel = %f\n",resultat_seq);
		
	}


	/* ******************************* */
	//		PRODUIT MATRICE VECTEUR 	/
	/* ******************************* */
	double* PMV_par;
	double* PMV_seq = NULL;

	//Tous les processus ont besoin du vecteur complet pour réaliser le produit matrice vecteur
	if(rank != 0){
		vecteur1_complet = (double *) malloc (nb_ligne*sizeof(double));
	}
	MPI_Bcast (vecteur1_complet, nb_ligne, MPI_DOUBLE, 0, MPI_COMM_WORLD);




	if(rank == 0){
		PMV_par = (double*) malloc (nb_ligne*sizeof(double));
		PMV_seq = (double*) malloc (nb_ligne*sizeof(double));
	}

	t1 = MPI_Wtime();
	for (i = 0; i < 1000; ++i)
	{
		produit_matrice_vecteur_parallele(matrice, vecteur1_complet, nb_ligne, nbligne_process, nb_col, PMV_par);
	}
	t2 = MPI_Wtime();

	if(rank == 0){
		printf("\n --- Produit Matrice Vecteur : matrice * vecteur1 ---\n");
		printf("\nTemps PMV parallèle : %f s\n", t2-t1);
		t1 = MPI_Wtime();
		for (i = 0; i < 1000; ++i)
		{
			produit_matrice_vecteur (matrice_complete, vecteur1_complet, nb_ligne, nb_col, PMV_seq);
		}
		t2 = MPI_Wtime();
		printf("\tResultat PMV parallèle =\n");
		//affiche(PMV_par, nb_ligne, nb_ligne);
		printf("Temps PMV séquentiel : %f s\n", t2-t1);
		printf("\tResultat PMV séquentiel = \n");
		//affiche(PMV_seq, nb_ligne, nb_ligne);
	}

	MPI_Barrier(MPI_COMM_WORLD);


	free(vecteur1);
	free(vecteur2);
	free(matrice);
	free(vecteur1_complet);
	if(rank == 0){
		free(vecteur2_complet);
		free(matrice_complete);
		free(PMV_seq);
		free(PMV_par);
	}

	MPI_Finalize();

	exit(0);
}