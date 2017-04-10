
	printf("\n");
	affiche(A, taille_matrice * taille_matrice, taille_matrice);
	printf("\n");





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




main() {

	printf("\n");
	int i;
	for(i = 0; i < taille_matrice; i++){
		if(valeurs_pro_i[i] == 0.0){
			//Affiches résultats valeurs propres réelles
			printf("valeur propre #%d = %g + %gi\n", i, valeurs_pro_r[i], valeurs_pro_i[i]);
			printf("\tvecteur propre droit associés =\n");
			affiche(vect_droit+i*taille_matrice, taille_matrice, taille_matrice);
			printf("\tvecteur propre gauche associés =\n");

			double* vect_droit_t =  transposee(vect_droit, taille_matrice);

			//Vérification A * vect = lambda * vect
			affiche(vect_gauche+i*taille_matrice, taille_matrice, taille_matrice);
			double* tmp_pmv = produit_matrice_vecteur(A, vect_droit+i*taille_matrice, taille_matrice, taille_matrice);
			double* tmp_mul = mul_vecteur_scalaire (valeurs_pro_r[i], vect_droit+i*taille_matrice, taille_matrice);
			printf("A*vect_pro = \n");
			affiche(tmp_pmv, taille_matrice, taille_matrice);
			printf("lambda * vect_pro = \n");
			affiche(tmp_mul, taille_matrice, taille_matrice);
		}
		printf("\n");
	}

	printf("\n--- Vecteurs propres droits (ROW-MAJOR) ---\n");
	affiche(vect_droit, taille_matrice * taille_matrice, taille_matrice);

	printf("\n--- Vecteurs propres gauches (ROW-MAJOR) ---\n");
	affiche(vect_gauche, taille_matrice * taille_matrice, taille_matrice);

	
}