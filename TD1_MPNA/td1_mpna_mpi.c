
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <mpi.h>
#include <sys/time.h> 


// struct decoupage_mpi {
// 	int debut;
// 	int fin;
// };

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


double produit_scalaire_parralle (double* vecteur1, double* vecteur2, int taille_vecteur, int debut, int fin, int rank){
	double produit_scalaire;
	int i;
	produit_scalaire = 0;

	printf("Process %d bien rentré dans fonction \n", rank);

	for(i=debut; i<fin; i++){
		produit_scalaire += vecteur1[i] * vecteur2[i];
		printf("Process %d : %d\n", rank, i);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	return produit_scalaire;
}


double* produit_matrice_vecteur_parrallele (double* matrice, double* vecteur, int nb_ligne, int nb_col, int debut, int fin){
	double* resultat;
	int i;

	resultat = (double*) malloc(nb_ligne*sizeof(double));

	for(i=debut; i<fin; i++){
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

	int size, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	nb_ligne = atoi(argv[1]);
	nb_col = atoi(argv[2]);


	if(rank == 0){
		printf("--- Debut du calcul des decoupage --- \n");
	}
	MPI_Barrier(MPI_COMM_WORLD);
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
	printf("\tRang %d :  nombre de lignes gérés %d\n",rank, nbligne_process);
	
	vecteur1 = (double*) malloc(nbligne_process*sizeof(double));
	vecteur2 = (double*) malloc(nbligne_process*sizeof(double));
	matrice  = (double*) malloc(nbligne_process*nb_col*sizeof(double));
	// printf("Nombres de mega-octets alloués %d \n", nb_ligne*nb_col*sizeof(double)/1000000);

	//Initialisation
	for(i=0; i<nbligne_process; i++){
		vecteur1[i] = 2; //(double) (rand()%100);
		vecteur2[i] = 2; //(double) (rand()%100);
		for(j=0; j<nb_col; j++){
			matrice[j + i*nb_col] = (double) i;
		}
	}

	// printf("Avant Send\n");
	// if(rank == 0){
	// 	for(i=0; i<size; i++){
	// 		MPI_Send(&vecteur1, 1, MPI_DOUBLE, i, 1000, MPI_COMM_WORLD);
	// 	}
	// }
	// else{
	// 	MPI_Recv(&vecteur1, 1, MPI_DOUBLE, 0, 1000, MPI_COMM_WORLD, NULL);
	// }

	// for(i=0; i<size; i++){
	// 	if(i==rank){
	// 		printf("Rang : %d, debut : %d, fin : %d \n",rank,debut,fin);
	// 		printf("vecteur1\n");
	// 		affiche(vecteur1, nb_ligne, nb_ligne);
	// 		printf("vecteur2\n");
	// 		affiche(vecteur2, nb_ligne, nb_ligne);		
	// 	}
	// 	MPI_Barrier(MPI_COMM_WORLD);
	// }
	// MPI_Barrier(MPI_COMM_WORLD);



	// struct decoupage_mpi decoupage = (struct decoupage_mpi) malloc (sizeof(struct decoupage_mpi));



	//---------------------------------------------------------------------------------------------//

	// En parallele

	// gettimeofday(&debut_calcul, NULL);
	// for (int i = 0; i < 100; ++i)
	// {
	// 	produit_scalaire_parralle(vecteur1, vecteur2, nb_ligne, debut, fin);
	// }
	// gettimeofday(&fin_calcul, NULL);
	// timersub(&fin_calcul, &debut_calcul, &duree_calcul);
	// printf("Temps produit scalaire parallèle : %f s\n", (double) (duree_calcul.tv_sec) + (duree_calcul.tv_usec / 1000000.0));


	// if(rank ==0){
	// 	printf("\nProduit scalaire vecteur1 vecteur2 : \n");
	// }
	// int produit_scalaire = produit_scalaire_parralle(vecteur1, vecteur2, nb_ligne, debut, fin, rank);
	// if(rank == 0){
	// 	printf("resultat %d\n",produit_scalaire);
	// }



	// gettimeofday(&debut_calcul, NULL);
	// for (int i = 0; i < 100; ++i)
	// {
	// 	produit_matrice_vecteur_parrallele(matrice, vecteur1, nb_ligne, nb_col, debut, fin);
	// }
	// gettimeofday(&fin_calcul, NULL);
	// timersub(&fin_calcul, &debut_calcul, &duree_calcul);

	// printf("Temps matrice vecteur parallèle : %f s\n", (double) (duree_calcul.tv_sec) + (duree_calcul.tv_usec / 1000000.0));

	free(vecteur1);
	free(vecteur2);
	free(matrice);

	MPI_Finalize();

	exit(0);
}