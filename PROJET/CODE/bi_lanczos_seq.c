
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
#include <lapacke.h>
// #include <cblas.h>

//Dans ce programme nos vecteurs sont stockés en ligne.


//Fonction qui calcule les valeurs/vecteurs propres d'un matrice de taille n
void calcul_valeurs_propres (double* matrice, int n,
							 double* vpropres_r, double* vpropres_i,
							  double* vect_gauche, double* vect_droit){

    int info;
    
    //Calcul des valeurs propres
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
void bi_lanczos(double* A, int m, int n, double* t, double *v, double* w){

	double alpha;
	double beta = 0;
	double delta = 0;

	// FILE* log = fopen("log.txt","a");

	//Vecteur initiaux
	v = (double *) malloc ((m+2)*n*sizeof(double));
	w = (double *) malloc ((m+2)*n*sizeof(double));
	t = (double *) malloc (m*m*sizeof(double));
	double* produit_A_v;
	double* produit_At_w;

	//Transposée de A
	double* At = transposee(A,n);

	//Initialisation des vecteurs initiaux
	int i;
	for(i=0; i<n; i++){
		v[i] = 0; //v_0
		w[i] = 0; //w_0
		v[i+n] = 0; //v_1
		w[i+n] = 0; //w_1
	}

	//pour que le (v_1,w_1) = 1
	v[n] = 1; 
	w[n] = 1;

	//Initialisation de la matrice T
	for(i=0; i<m*m; i++){
		t[i] = 0;
	}

	//Check si les vecteurs initiaux respectent la propriété (v_1,w_1) = 1
	if(produit_scalaire(v+n, w+n, n) != 1.){
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
	// affiche(transposee(t, m), m*m, m);
}

//Fonction de test des différentes fonctions implémentées
void test(){

	int n = 5;
	int m = 8;
	double* vecteur1 = (double *) malloc (n*sizeof(double));
	double* vecteur2 = (double *) malloc (n*sizeof(double));
	double* vecteur3 = (double *) malloc (n*sizeof(double));
	int i;
	for(i=0; i<n; i++){
	       vecteur1[i] = i;
	       vecteur2[i] = 2*i;
	       vecteur3[i] = 3*i;
	}


	// printf("Vecteur1 :\n");
	// affiche(vecteur1, n, n);
	// printf("Vecteur2 :\n");
	// affiche(vecteur2, n, n);
	// printf("Vecteur3 :\n");
	// affiche(vecteur3, n, n);

	// //Test vecteur * scalaire
	// double* result = mul_vecteur_scalaire(4, vecteur1, n);
	// printf("--- TEST vect*scal ---\n  Vecteur1 * 4 =\n");
	// affiche(result, n, n);

	// //Test vecteur - vecteur
	// diff_vecteurs(vecteur1, vecteur1, result, n);
	// printf("--- TEST vect-vect ---\n  Vecteur1 - Vecteur1 =\n");
	// affiche(result, n, n);

	// //Test vecteur / scalaire
	// div_vecteur_scalaire(vecteur3, 3, n);

	// printf("--- TEST vect/scal ---\n  Vecteur3 / 3 =\n");
	// affiche(vecteur3, n, n);

	double * matrice1 = (double *) malloc (n*n*sizeof(double));
	double * matrice2 = (double *) malloc (n*n*sizeof(double));

	int j;
	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			if(j==i){
				matrice1[i + j*n] = rand()%3 + 0.;
			}
			else{
				matrice1[i + j*n] = 0;
			}
			matrice2[i + j*n] = 2;
		}
	}
	// printf("matrice1 : \n");
	// affiche(matrice1, n*n, n);
	// printf("matrice2 : \n");
	// affiche(matrice2, n*n, n);
	// double* resultat = produit_matrice_matrice(matrice1, matrice2, n, n, n, n);
	// printf("resultat : \n");
	// affiche(resultat, n*n, n);


	//Tests transposée générale avec matrice (n*m)
	double * matrice3 = (double *) malloc (n*m*sizeof(double));
	for(i=0; i < n*m; i++){
		matrice3[i] = rand()%9 + 0.;
	}
	printf("matrice3 : \n");
	affiche(matrice3, n*m, m);
	double* matrice3_t = transposee2(matrice3, n, m);
	printf("transposée matrice3 :\n");
	affiche(matrice3_t, n*m, n);

	free(vecteur1);
	free(vecteur2);
	free(vecteur3);
	free(matrice1);
	free(matrice2);
	free(matrice3);
	// free(resultat);
}


int main(int argc, char const *argv[]){

	if(argc < 3){
		fprintf(stderr, "Erreur arguments:\n\tUsage : %s <taille_matrice> <taille_ss_espace_krylov>\n", argv[0]);
		exit(-1);
	}

	int taille_ss_espace_krylov = atoi(argv[2]);
	int taille_matrice = atoi(argv[1]);

	fprintf(stderr, "taille_matrice = %d, taille_ss_espace_krylov = %d\n",taille_matrice, taille_ss_espace_krylov);

	if (taille_ss_espace_krylov > taille_matrice){
		fprintf(stderr, "Erreur arguments:\n\tUsage : La <taille_ss_espace_krylov>  doit être inférieur ou égale à la <taille_matrice>\n");
		exit(-1);
	}



	struct timeval debut_calcul, fin_calcul, duree_calcul;

	//Allocation et initialisation de la matrice à projeter
	double* A = (double *) malloc (taille_matrice*taille_matrice*sizeof(double));
	double* v;
	double* w;
	double* t;

	int i;
	for (i = 0; i < taille_matrice*taille_matrice; i++)
	{
		A[i] = rand()%80 + 10.;
	}

	FILE* output = fopen("OUTPUT/temps_execution.txt","a+");

	if(output == NULL){
		fprintf(stderr, "Erreur d'ouverture du fichier OUTPUT/temps_execution.txt\n");
		exit(-1);
	}

	if(taille_matrice <= 10){
		affiche(A, taille_matrice*taille_matrice, taille_matrice);
		printf("transposee : \n");
		affiche(transposee(A,taille_matrice), taille_matrice*taille_matrice, taille_matrice);
	}


	gettimeofday(&debut_calcul, NULL);
	bi_lanczos(A,taille_ss_espace_krylov, taille_matrice, t, v, w);
	gettimeofday(&fin_calcul, NULL);
	timersub(&fin_calcul, &debut_calcul, &duree_calcul);
	fprintf(output, "%d %f\n", taille_matrice,
			(double) (duree_calcul.tv_sec) + (duree_calcul.tv_usec / 1000000.0));

	return 0;
}