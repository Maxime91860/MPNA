
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

double* mul_vecteur_constante (double constante, double* vecteur, int n){
	double* resultat;
	int i;

	resultat = (double*) malloc(n*sizeof(double));

	for(i=0; i<n; i++){
		resultat[i] = vecteur[i] * constante;
	}

	return resultat;
}

void diff_vecteurs (double* vecteur1, double* vecteur2, double* resultat, int n){
	int i;

	for(i=0; i<n; i++){
		resultat[i] = vecteur1[i] - vecteur2[i];
	}
}

void vecteur_div_constante(double* vecteur, double constante, int n){

	int i;
	for (i = 0; i < n; i++)
	{
		vecteur[i] = vecteur[i] / constante;
	}
}

//Transposée d'une matrice carrée
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

void bi_lanczos(double* A, int m, int n){

	double alpha;
	double beta = 0;
	double delta = 0;

	FILE* log = fopen("log.txt","a");


	//Vecteur initiaux
	double* v = (double *) malloc ((m+2)*n*sizeof(double));
	double* w = (double *) malloc ((m+2)* n*sizeof(double));
	double* t = (double *) malloc (m*m*sizeof(double));
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

	int j;
	for(j=1; j <= m; j++){
		produit_A_v = produit_matrice_vecteur (A, v+(j*n), n, n);
		produit_At_w = produit_matrice_vecteur (At, w+(j*n), n, n);

		alpha = produit_scalaire(produit_A_v, w+(j*n), n);

		fprintf(log, "ITERATION #%d : alpha = %f, beta = %f, delta = %f\n", j, alpha, beta, delta);

		double* res1;
		double* res2;
		double* res3 = (double *) malloc (n*sizeof(double));

		//Calcul de v_j+1
		res1 = mul_vecteur_constante(alpha, v+(j*n), n);
		res2 = mul_vecteur_constante(beta, v+((j-1)*n), n);
		diff_vecteurs(res1, res2, res3, n);
		diff_vecteurs(produit_A_v, res3, v+((j+1)*n), n);

		//Calcul de w_j+1
		res1 = mul_vecteur_constante(alpha, w+(j*n), n);
		res2 = mul_vecteur_constante(delta, w+((j-1)*n), n);
		diff_vecteurs(res1, res2, res3, n);
		diff_vecteurs(produit_At_w, res3, w+((j+1)*n), n);

		delta = sqrt(abs(produit_scalaire(v+((j+1)*n), w+((j+1)*n), n)));
		if(delta == 0){
			fprintf(stderr, "Erreur numérique du calcul de delta#%d.\n",j);
			exit(-1);
		}

		beta = produit_scalaire (v+((j+1)*n), w+((j+1)*n), n) / delta;


		vecteur_div_constante(v+((j+1)*n), beta, n);
		vecteur_div_constante(w+((j+1)*n), delta, n);

		//Remplissage de T
		printf("alpha = %g\n",alpha);
		t[j-1+(j-1)*m] = alpha;
		if(j != m){
			t[j-1+(j)*m] = beta;
			t[j+(j-1)*m] = delta;
		}
	}

	printf("V =\n");
	affiche(v, n*(m+1), n);
	printf("W =\n");
	affiche(w, n*(m+1), n);
	printf("T =\n");
	affiche(transposee(t, m), m*m, m);

}

void test(){
	int n = 10;
	double* vecteur1 = (double *) malloc (n*sizeof(double));
	double* vecteur2 = (double *) malloc (n*sizeof(double));
	double* vecteur3 = (double *) malloc (n*sizeof(double));

	int i;
	for(i=0; i<n; i++){
		vecteur1[i] = i;
		vecteur2[i] = 2*i;
		vecteur3[i] = 3*i;
	}


	printf("Vecteur1 :\n");
	affiche(vecteur1, n, n);
	printf("Vecteur2 :\n");
	affiche(vecteur2, n, n);
	printf("Vecteur3 :\n");
	affiche(vecteur3, n, n);

	//Test vecteur * scalaire
	double* result = mul_vecteur_constante(4, vecteur1, n);
	printf("--- TEST vect*scal ---\n  Vecteur1 * 4 =\n");
	affiche(result, n, n);

	//Test vecteur - vecteur
	diff_vecteurs(vecteur1, vecteur1, result, n);
	printf("--- TEST vect-vect ---\n  Vecteur1 - Vecteur1 =\n");
	affiche(result, n, n);

	//Test vecteur * scalaire
	vecteur_div_constante(vecteur3, 3, n);
	printf("--- TEST vect/scal ---\n  Vecteur3 / 3 =\n");
	affiche(vecteur3, n, n);

}

int main(int argc, char const *argv[])
{

	if(argc < 3){
		fprintf(stderr, "Erreur :\nUsage : %s <taille_matrice> <taille_ss_espace_krylov>\n", argv[0]);
		exit(-1);
	}

	
	int taille_ss_espace_krylov = atoi(argv[2]);
	int taille_matrice = atoi(argv[1]);
	double* A = (double *) malloc (taille_matrice*taille_matrice*sizeof(double));

	int i;
	for (i = 0; i < taille_matrice*taille_matrice; i++)
	{
		A[i] = rand()%80 + 10.;
	}

	if(taille_matrice <= 10){
		affiche(A, taille_matrice*taille_matrice, taille_matrice);
		printf("transposee : \n");
		affiche(transposee(A,taille_matrice), taille_matrice*taille_matrice, taille_matrice);
	}


	

	// test();

	bi_lanczos(A,taille_ss_espace_krylov, taille_matrice);

	return 0;
}