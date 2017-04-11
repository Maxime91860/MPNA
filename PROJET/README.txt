
Se placer dans le répertoire courant de ce README.


Pour compiler les programmes :

-> make 


Pour executer :

-> ./CODE/bi_lanczos_seq.pgr <fichier_matrix_market> <taille_ss_espace_krylov>


Pour compiler et executer un exemple (la matrice bp_400.mtx non-symétrique réelle) :

-> make exemple


Nettoyer les binaires :

-> make clean

La sortie du programme est la projection de la matrice donnée en entrée avec la méthode bi-lanczos,
ainsi que les valeurs propres de cette matrice, et les valeurs propres de sa projection.

Notre implementation de la méthode bi-lanczos est dans le fichier "bi_lanczos_seq.c" dans la fonction :
	void bi_lanczos(double* A, int m, int n, double* t, double *v, double* w)

L'appel à la bibliothèque LAPACK est dans la fonction :
	void calcul_valeurs_propres (double* matrice, int n, double* vpropres_r, double* vpropres_i, double* vect_gauche, double* vect_droit)

