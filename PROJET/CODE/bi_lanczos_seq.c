
/**********************************************************************************/
/* Implémentation séquentielle de l'algorithme de bi-orthogonalisation de Lanczos */
/*				 Amira AKLOUL - Maxime KERMARQUER --- M2 CHPS					  */
/**********************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/time.h> 
#include <math.h>
#include <string.h>
#include <sys/time.h> 
#include "mmio.h" //Fonctions qui parsent les fichiers de matrix market
#include <lapacke.h>
// #include <cblas.h>

//Dans ce programme nos vecteurs sont stockés en ligne.


//Fonction qui lit un fichier .mtx et renvoie la matrice du fichier
double* lis_matrice (char* file_mtx_name, int* taille_matrice){

    MM_typecode matcode;
    FILE *file;
    int nb_ligne, nb_col, non_zero;   
    int i;
	int *II;
    int *J;
    double *val;


    fprintf(stderr, "-- Debut lecture matrice %s --\n",file_mtx_name);
    if ((file = fopen(file_mtx_name, "r")) == NULL) 
        exit(1);

    if (mm_read_banner(file, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }

    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */
    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
            mm_is_sparse(matcode) )
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    //Récupération des dimensions de la matrice et du nombres d'éléments non-nuls
    if (mm_read_mtx_crd_size(file, &nb_ligne, &nb_col, &non_zero) !=0)
        exit(1);

    II = (int *) malloc(non_zero * sizeof(int));
    J = (int *) malloc(non_zero * sizeof(int));
    val = (double *) malloc(non_zero * sizeof(double));

    for (i=0; i<non_zero; i++)
    {
        fscanf(file, "%d %d %lg\n", &II[i], &J[i], &val[i]);
        II[i]--;  /* adjust from 1-based to 0-based */
        J[i]--;
    }

    fclose(file);

    //Passage de stockage creux à stockage plein
    double* matrice = (double*) malloc(nb_ligne * nb_col * sizeof(double));
    memset(matrice, 0, nb_ligne * nb_col * sizeof(double));
    for (i=0; i<non_zero; i++)
    {
    	matrice [II[i]*nb_ligne + J[i]] = val[i];
    }   

    if(nb_ligne != nb_col) {
    	fprintf(stderr, "Erreur : la matrice n'est pas carrée.\n");
    	exit(-1);
    }

    fprintf(stderr, "-- Fin lecture matrice %s --\n",file_mtx_name);
    fprintf(stderr, "La matrice %s est de dimension %dx%d et a %d élements non-nuls.\n",file_mtx_name, nb_ligne, nb_col, non_zero);

    *taille_matrice = nb_ligne;
    return matrice;
}

//Fonction qui renvoie une matrice aléatoire
double* init_matrice (int taille_matrice){
	int i;
	double* A = (double *) malloc (taille_matrice*taille_matrice*sizeof(double));
	for (i = 0; i < taille_matrice*taille_matrice; i++)
	{
		A[i] = rand()%80 + 10.;
	}
	return A;
}

//Fonction qui calcule les valeurs/vecteurs propres d'un matrice de taille n
void calcul_valeurs_propres (double* matrice, int n, double* vpropres_r, double* vpropres_i, double* vect_gauche, double* vect_droit){

    int info;
    
    //Calcul des valeurs propres
    //Details des arguments :
    //#1 Comment la matrice est stocker
    //#2 Calcul ou non des vecteurs propres à gauche
    //#3 Calcul ou non des vecteurs propres à droite
    //#4 La dimension de la matrice
    //#5 La matrice
    //#6 La dimension de la matrice
    //#7 Les parties reelles des valeurs propres
    //#8 Les parties imaginaires des valeurs propres
    //#9 Les vecteurs propres à gauche
    //#10 La taille des vecteurs propres à gauche
    //#11 Les vecteurs propres à droitre
    //#12 La taille des vecteurs propres à droite
    info = LAPACKE_dgeev( LAPACK_ROW_MAJOR, 'V', 'V', n, matrice, n, vpropres_r, vpropres_i, vect_gauche, n, vect_droit, n );

    //Vérification du succès du calcul
    if( info > 0 ) {
            printf( "Erreur : Calcul valeurs propres.\n" );
            exit( 1 );
    }
}

//Fonction d'affichage d'un vecteur ou d'une matrice
void affiche(double* tab, int N, int width){

	int i;
	printf("(");
	for(i=0; i<N; i++){
		if((i+1)%width == 0){
			printf("%g)\n",tab[i]);
			if(i != (N-1)){
				printf("(");
			}
		}
		else{
			printf("%g ",tab[i]);
		}
	}
}

//Calcul de la norme 2 d'un vecteur
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

//Calcul du produit scalaire de deux vecteurs
double produit_scalaire (double* vecteur1, double* vecteur2, int taille_vecteur){

	double produit_scalaire;
	int i;
	produit_scalaire = 0;

	for(i=0; i<taille_vecteur; i++){
		produit_scalaire += vecteur1[i] * vecteur2[i];
	}

	return produit_scalaire;
}

//Calcul du produit matrice par vecteur
double* produit_matrice_vecteur (double* matrice, double* vecteur, int nb_ligne, int nb_col){

	double* resultat;
	int i;

	resultat = (double*) malloc(nb_ligne*sizeof(double));

	for(i=0; i<nb_ligne; i++){
		resultat[i] = produit_scalaire(matrice+(i*nb_col), vecteur, nb_col);
	}

	return resultat;
}

//Calcul du produit d'un vecteur par un scalaire
double* mul_vecteur_scalaire (double scalaire, double* vecteur, int n){

	double* resultat;
	int i;

	resultat = (double*) malloc(n*sizeof(double));

	for(i=0; i<n; i++){
		resultat[i] = vecteur[i] * scalaire;
	}

	return resultat;
}

//Calcul de la division d'un vecteur par un scalaire
void div_vecteur_scalaire (double* vecteur, double scalaire, int n){

	int i;
	for (i = 0; i < n; i++)
	{
		vecteur[i] = vecteur[i] / scalaire;
	}
}

//Calcul de la différence entre deux vecteurs
void diff_vecteurs (double* vecteur1, double* vecteur2, double* resultat, int n){

	int i;
	for(i=0; i<n; i++){
		resultat[i] = vecteur1[i] - vecteur2[i];
	}
}

//Calcul de la transposée d'une matrice carrée
double* transposee (double* A, int n){

	double * At = (double *) malloc (n*n*sizeof(double));
	int i, j; 
	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			At[j+i*n] = A[i+j*n];
		}
	}
	return At;
}

//Calcul de la transposée d'une matrice (n*m)
double* transposee2 (double* A, int n, int m){

	double * At = (double *) malloc (m*n*sizeof(double));
	int i, j; 
	for(i=0; i<m; i++){
		for(j=0; j<n; j++){
			At[j+i*n] = A[i+j*m];
		}
	}
	return At;
}

//Calcul du produit matriciel de la matrice1 (n1*m1) par matrice2 (n2*m2) (à valider)
double* produit_matrice_matrice (double* matrice1, double* matrice2, int n1, int m1, int n2, int m2){

	if(m1 != n2){
		fprintf(stderr, "Erreur fonctions produit_matrice_matrice:\n--> m1 != n2\n");
		return NULL;
	}

	int i, j;
	// double* tmat2 = transposee(matrice2, n2);
	double* resultat = (double *) malloc (n1*m2*sizeof(double));
	for(i=0; i<n1; i++){
		for(j=0; j<m2; j++){
			resultat [i+ j*n1] = produit_scalaire(matrice1+(i*m1), matrice2+(j*m2), n2);
		}
	}
	return resultat;
}

//Méthode de projection biorthogonale de Lanczos
//On va creer :
//	- deux bases V et W
//	- la projection de la matrice A, dans T
// V sera la base pour A
// W sera la base pour la transposée de A, A^t
// m dimension du sous-espace de krylov
// n dimension de la matrice
void bi_lanczos(double* A, int m, int n, double* t, double *v, double* w){

	double alpha;
	double beta = 0;
	double delta = 0;

	// FILE* log = fopen("log.txt","a");

	double* produit_A_v;
	double* produit_At_w;

	//Transposée de A
	double* At = transposee(A,n);

	//Initialisation des vecteurs initiaux
	int i;
	for(i=0; i<n; i++){
		v[i] = 0; //v_0
		w[i] = 0; //w_0
		v[i+n] = rand()%10; //v_1
		w[i+n] = rand()%10; //w_1
	}

	//Pour que le (v_1,w_1) = 1
	double ps = produit_scalaire(w+n, v+n, n);
	for(i = 0; i<n; i++){
		w[i+n] = w[i+n] / ps;
	}
	// printf("v_1\n");
	// affiche(v+n, n, n);
	// printf("w_1\n");
	// affiche(w+n, n, n);
	// printf("w_1.v_1 = %f\n", produit_scalaire(w+n, v+n, n));

	//Initialisation de la matrice T
	for(i=0; i<m*m; i++){
		t[i] = 0;
	}

	//Check si les vecteurs initiaux respectent la propriété (v_1,w_1) = 1
	if(fabs(produit_scalaire(v+n, w+n, n) - 1) > 1.e-5){
		fprintf(stderr, "Le produit scalaire des vecteurs initiaux doit être égal à 1.\n");
		exit(-1);
	}


	/*---------------------------------------------------*/
	/*--- PROCEDURE DE BIORTHOGONALISATION DE LANCZOS ---*/
	/*---------------------------------------------------*/ 
	int j;
	for(j=1; j <= m; j++){
		//Calcul des produits matrice * vecteur
		produit_A_v = produit_matrice_vecteur (A, v+(j*n), n, n);
		produit_At_w = produit_matrice_vecteur (At, w+(j*n), n, n);

		alpha = produit_scalaire(produit_A_v, w+(j*n), n);

		// fprintf(log, "ITERATION #%d : alpha = %f, beta = %f, delta = %f\n", j, alpha, beta, delta);

		double* res1;
		double* res2;
		double* res3 = (double *) malloc (n*sizeof(double));

		//Calcul de v_(j+1)
		res1 = mul_vecteur_scalaire(alpha, v+(j*n), n);
		res2 = mul_vecteur_scalaire(beta, v+((j-1)*n), n);
		diff_vecteurs(res1, res2, res3, n);
		diff_vecteurs(produit_A_v, res3, v+((j+1)*n), n);

		//Calcul de w_(j+1)
		res1 = mul_vecteur_scalaire(alpha, w+(j*n), n);
		res2 = mul_vecteur_scalaire(delta, w+((j-1)*n), n);
		diff_vecteurs(res1, res2, res3, n);
		diff_vecteurs(produit_At_w, res3, w+((j+1)*n), n);

		delta = sqrt(abs(produit_scalaire(v+((j+1)*n), w+((j+1)*n), n)));
		if(delta == 0){
			fprintf(stderr, "Erreur numérique du calcul de delta à l'itération #%d.\n",j);
			exit(-1);
		}

		beta = produit_scalaire (v+((j+1)*n), w+((j+1)*n), n) / delta;

		//MaJ de w_(j+1)
		div_vecteur_scalaire (w+((j+1)*n), beta, n);

		//MaJ de v_(j+1)
		div_vecteur_scalaire (v+((j+1)*n), delta, n);

		//Remplissage de T
		t[j-1+(j-1)*m] = alpha;
		if(j != m){
			t[j-1+(j)*m] = beta; 
			t[j+(j-1)*m] = delta;
		}
	}



	// printf("V =\n");
	// affiche(v, n*(m+1), n);
	// printf("W =\n");
	// affiche(w, n*(m+1), n);
	// printf("T =\n");
	// affiche(t, m*m, m);
}


int main(int argc, char* argv[]){

	//Check arguments
	if(argc < 3){
		fprintf(stderr, "Erreur arguments:\n\tUsage : %s <fichier_matrix_market> <taille_ss_espace_krylov>\n", argv[0]);
		exit(-1);
	}

	int taille_ss_espace_krylov = atoi(argv[2]);
	int taille_matrice = 40;

	//Matrice à projeter
	double* A = lis_matrice(argv[1], &taille_matrice); 
	// double* A = init_matrice(40);

	//Check arguments 2
	if (taille_ss_espace_krylov > taille_matrice && taille_matrice > 5){
		fprintf(stderr, "Erreur arguments:\n\tUsage : La <taille_ss_espace_krylov> doit être inférieure ou égale à la <taille_matrice>.\n");
		exit(-1);
	}


	//v contiendra l'ensemble des vecteurs v, v[0] renvoie au vecteur v_0 par exemple
	double* v;

	//w contiendra l'ensemble des vecteurs w
	double* w;

	//t est la projection de la matrice A, t est de dimension m = taille_ss_espace_krylov
	double* t;


	struct timeval debut_calcul, fin_calcul, duree_calcul;
	FILE* output = fopen("OUTPUT/temps_execution.txt","a+");
	if(output == NULL){
		fprintf(stderr, "Erreur d'ouverture du fichier OUTPUT/temps_execution.txt\n");
		exit(-1);
	}

	if(taille_matrice <= 10){
		printf("\n--- Matrice A à projeter ---\n");
		affiche(A, taille_matrice*taille_matrice, taille_matrice);
		// printf("\n--- Transposée A^t ---\n");
		// affiche(transposee(A,taille_matrice), taille_matrice*taille_matrice, taille_matrice);
	}


	//Allocations vecteur initiaux
	v = (double *) malloc ((taille_ss_espace_krylov+2)*taille_matrice*sizeof(double));
	w = (double *) malloc ((taille_ss_espace_krylov+2)*taille_matrice*sizeof(double));

	//Allocation de la projection de matrice A
	t = (double *) malloc (taille_ss_espace_krylov*taille_ss_espace_krylov*sizeof(double));


	gettimeofday(&debut_calcul, NULL);
	if(taille_ss_espace_krylov > 5)
		bi_lanczos(A,taille_ss_espace_krylov, taille_matrice, t, v, w);
	gettimeofday(&fin_calcul, NULL);

	if (taille_ss_espace_krylov <= 10 && taille_ss_espace_krylov < taille_matrice){
		printf("\n--- Projection T de la matrice A ---\n");
		affiche(t, taille_ss_espace_krylov*taille_ss_espace_krylov, taille_ss_espace_krylov);
	}

	//Mesure du temps d'execution de la projection, et ecriture dans un fichier de sortie utile pour les plots
	timersub(&fin_calcul, &debut_calcul, &duree_calcul);
	fprintf(output, "%d %f\n", taille_matrice,
			(double) (duree_calcul.tv_sec) + (duree_calcul.tv_usec / 1000000.0));

	//Calcul des valeurs/vecteurs propres exacts
	double *vect_droit = (double *) malloc (taille_matrice*taille_matrice*sizeof(double));
	double *vect_gauche = (double *) malloc (taille_matrice*taille_matrice*sizeof(double));;
	double *valeurs_pro_r = (double *) malloc (taille_matrice*sizeof(double));
	double *valeurs_pro_i = (double *) malloc (taille_matrice*sizeof(double));;
	
	calcul_valeurs_propres(A, taille_matrice, valeurs_pro_r, valeurs_pro_i, vect_droit, vect_gauche);


	//Calcul des valeurs/vecteurs propres approchés
	double *vect_droit_app = (double *) malloc (taille_ss_espace_krylov*taille_ss_espace_krylov*sizeof(double));
	double *vect_gauche_app = (double *) malloc (taille_ss_espace_krylov*taille_ss_espace_krylov*sizeof(double));;
	double *valeurs_pro_r_app = (double *) malloc (taille_ss_espace_krylov*sizeof(double));
	double *valeurs_pro_i_app = (double *) malloc (taille_ss_espace_krylov*sizeof(double));;
	
	calcul_valeurs_propres(t, taille_ss_espace_krylov, valeurs_pro_r_app, valeurs_pro_i_app, vect_droit_app, vect_gauche_app);

	
	//v+taille_matrice pour ne pas prendre en compte le v_0
	// vect_droit_app = produit_matrice_vecteur(v+taille_matrice, vect_droit_app, taille_matrice, taille_ss_espace_krylov);
	// vect_gauche_app = produit_matrice_vecteur(w+taille_matrice, vect_gauche_app, taille_matrice, taille_ss_espace_krylov);



	if (taille_matrice > taille_ss_espace_krylov)
		taille_matrice = taille_ss_espace_krylov;



	printf("\n----------------------------------------------------\n--- CALCULS DES VALEURS PROPRES SUR LA MATRICE A ---\n----------------------------------------------------\n");
	printf("\n---  Valeurs propres parties réelles ---\n");
	affiche(valeurs_pro_r, taille_matrice, taille_matrice);

	printf("\n---  Valeurs propres parties imaginaires ---\n");
	affiche(valeurs_pro_i, taille_matrice, taille_matrice);

	// printf("\n--- Vecteurs propres droits (ROW-MAJOR) ---\n");
	// affiche(vect_droit, taille_matrice * taille_matrice, taille_matrice);

	// printf("\n--- Vecteurs propres gauches (ROW-MAJOR) ---\n");
	// affiche(vect_gauche, taille_matrice * taille_matrice, taille_matrice);


	//Le if pour le cas où on n'effectue pas la projection.
	if(taille_ss_espace_krylov <= taille_matrice){
		printf("\n-------------------------------------------------------\n--- CALCULS DES VALEURS PROPRES SUR LA PROJECTION T ---\n-------------------------------------------------------\n");
		printf("\n---  Valeurs propres parties réelles ---\n");
		affiche(valeurs_pro_r_app, taille_ss_espace_krylov, taille_ss_espace_krylov);

		printf("\n---  Valeurs propres parties imaginaires ---\n");
		affiche(valeurs_pro_i_app, taille_ss_espace_krylov, taille_ss_espace_krylov);

		// printf("\n--- Vecteurs propres droits (ROW-MAJOR) ---\n");
		// affiche(vect_droit_app, taille_matrice * taille_matrice, taille_matrice);

		// printf("\n--- Vecteurs propres gauches (ROW-MAJOR) ---\n");
		// affiche(vect_gauche_app, taille_matrice * taille_matrice, taille_matrice);
	}
	return 0;
}