#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include <hdf5.h>
#include <hdf5_hl.h>

#define VERBOSE
//#define RANK 2
#define VAL2D(arr, xx, yy) (arr[(xx)+width*(yy)])

void ecrire(double* cur, int width, int height, const char* file, const char* set)
{
	hsize_t dim[2];
	dim[0]=width;
	dim[1]=height;
	hid_t file_id=H5Fopen(file, H5F_ACC_RDWR, H5P_DEFAULT);

	H5LTmake_dataset_double(file_id, set, 2, dim, cur);

	H5Fclose(file_id);
}

void lire(double* cur, int width, int height, const char* file, const char* set)

{
	hid_t file_id=H5Fopen(file, H5F_ACC_RDWR, H5P_DEFAULT);
	
	if(file_id<0)
	{
		fprintf(stderr, "Unreachable HDF5 file : %s\n", file);
		exit(EXIT_FAILURE);
	}
	
	H5LTread_dataset_double(file_id, set, cur);
	
	H5Fclose(file_id);
}

// Calcul des moyennes

void moyennes(double* cur, int width, int height, int n, int rank, const char* file)
{
	int xx, yy, Q, R, ideb, ifin;
	int i, j;
	int tag = 1000;
	MPI_Status sta;
	double moyenne_lignes = 0, moyenne_colonnes = 0, moyenne_globale;
	
	Q = height / n;
	R = height % n;
	
		if(rank < R){
			ideb = rank * (1+Q);
			ifin = ideb + Q;
		}
		else{
		ideb = R * (1+Q) + ((rank - R) * Q);
		ifin = ideb + Q -1;	
		}
	
	double *moy_x = (double *)malloc((1+ifin-ideb)*sizeof(double));
	double somme;
	
	for(i = ideb; i<=ifin; i++){
		somme = 0;
		for(xx=0; xx<width; xx++){
			somme = somme + VAL2D(cur, xx, i);
		}
		
		moy_x[i-ideb] = somme/width;

	}

		if (rank != 0){
			for (i=0; i <= (ifin-ideb); i++){
				MPI_Send(&moy_x[i], 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
				}
		}
		else if(rank==0){
			double *moyenne_x;
			moyenne_x = (double *)malloc(sizeof(double) * height);
				for (j=0; j<=(ifin-ideb); j++){
					moyenne_x[j] = moy_x[j];
				}
			
				for(int t = 1 ; t < n ; t++ ) {
					
					if (R>1){
					if (t < R){
						for(i = 0 ; i <= (ifin-ideb) ; i++ ){
						int index = (((ifin+1)*t)+i);
						MPI_Recv(&moyenne_x[index], 1, MPI_DOUBLE, t, tag, MPI_COMM_WORLD, &sta);
						}
					}
					else{
						for(i = 0 ; i <= ((ifin-ideb)-1) ; i++ ){
						int index= (((ifin+1)*R)+((t-R)*(ifin-ideb))+i);
						MPI_Recv(&moyenne_x[index], 1, MPI_DOUBLE, t, tag, MPI_COMM_WORLD, &sta);
						}
					}
					}
					else if (R==1){
						for(i = 0 ; i <= ((ifin-ideb)-1) ; i++ ){
						int index= (ifin+1+((t-1)*(ifin-ideb))+i);
						MPI_Recv(&moyenne_x[index], 1, MPI_DOUBLE, t, tag, MPI_COMM_WORLD, &sta);
						}
						}
					else{
						for(i = 0 ; i <= (ifin-ideb) ; i++ ){
						MPI_Recv(&moyenne_x[(((ifin+1)*t)+i)], 1, MPI_DOUBLE, t, tag, MPI_COMM_WORLD, &sta);
						}
					}
				}
				
				for(i=0; i<height; i++){
					#ifdef VERBOSE
					printf("La moyenne de la ligne %d = %f\n",i, moyenne_x[i]);
					#endif  
					moyenne_lignes = moyenne_lignes + moyenne_x[i];
				}
				moyenne_lignes = moyenne_lignes / height;				
				free(moyenne_x);
			}
			
				if (rank==0){
				
				//printf("La moyenne des lignes = %f\n", moyenne_lignes);
				printf("La moyenne des lignes = %f\n", moyenne_lignes);
				ecrire(&moyenne_lignes, 1, height, file, "/moyennesLignes");
				
				}
		
		// Calcule des moyennes des colonnes
		
	Q = width / n;
	R = width % n;
	
		if(rank < R){
			ideb = rank * (1+Q);
			ifin = ideb + Q;
		}
		else{
		ideb = R * (1+Q) + ((rank - R) * Q);
		ifin = ideb + Q -1;	
		}
	
	double *moy_y = (double *)malloc((1+ifin-ideb)*sizeof(double));	
	
	for(i = ideb; i<=ifin; i++){
		somme = 0;
		for(xx=0; xx<height; xx++){
			somme = somme + VAL2D(cur, i, xx);
		}
		
		moy_y[i-ideb] = somme/height;
	}

		if (rank != 0){
			for (i=0; i <= (ifin-ideb); i++){
				MPI_Send(&moy_y[i], 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
				}
		}
		else if(rank==0){
			double *moyenne_y;
			moyenne_y = (double *)malloc(sizeof(double) * height);
				for (j=0; j<=(ifin-ideb); j++){
					moyenne_y[j] = moy_y[j];
				}
			
				for(int t = 1 ; t < n ; t++ ) {
					
					if (R>1){
					if (t < R){
						for(i = 0 ; i <= (ifin-ideb) ; i++ ){
						int index = (((ifin+1)*t)+i);
						MPI_Recv(&moyenne_y[index], 1, MPI_DOUBLE, t, tag, MPI_COMM_WORLD, &sta);
						}
					}
					else{
						for(i = 0 ; i <= ((ifin-ideb)-1) ; i++ ){
						int index= (((ifin+1)*R)+((t-R)*(ifin-ideb))+i);
						MPI_Recv(&moyenne_y[index], 1, MPI_DOUBLE, t, tag, MPI_COMM_WORLD, &sta);
						}
					}
					}
					else if (R==1){
						for(i = 0 ; i <= ((ifin-ideb)-1) ; i++ ){
						int index= (ifin+1+((t-1)*(ifin-ideb))+i);
						MPI_Recv(&moyenne_y[index], 1, MPI_DOUBLE, t, tag, MPI_COMM_WORLD, &sta);
						}
						}
					else{
						for(i = 0 ; i <= (ifin-ideb) ; i++ ){
						MPI_Recv(&moyenne_y[(((ifin+1)*t)+i)], 1, MPI_DOUBLE, t, tag, MPI_COMM_WORLD, &sta);
						}
					}
				}
				
				for(i=0; i<width; i++){
					#ifdef VERBOSE
					printf("La moyenne de la colonne %d = %f\n",i, moyenne_y[i]);
					#endif  
					moyenne_colonnes = moyenne_colonnes + moyenne_y[i];
				}
				moyenne_colonnes = moyenne_colonnes / width;
				#ifdef VERBOSE
				printf("La moyenne des colonnes = %f\n", moyenne_colonnes);
				#endif
				ecrire(&moyenne_colonnes, width, 1, file, "/moyennesColones");
				free(moyenne_y);
		}
	
	if(rank==0){
	moyenne_globale = (moyenne_lignes + moyenne_colonnes) / 2;
	#ifdef VERBOSE
	printf("La moyenne globale = %f\n", moyenne_globale);
	#endif
	ecrire(&moyenne_globale, 1, 1, file, "/moyenneGlobale");
	}
	free(moy_x);
	free(moy_y);
}

void laplacien(double* cur, int width, int height, const char* file)
{
	int xx, yy;
    double * laplacien = (double*)malloc((width-2)*(height-2)*sizeof(double));
	
	for (int j = 0; j<(height-2); j++){
		for(int i=0; i<(width-2); i++){
			laplacien[i+(j*(width-2))] = (-4)* VAL2D(cur, j+1, i+1) + VAL2D(cur, j+2, i+1)
									  + VAL2D(cur, j, i+1) + VAL2D(cur, j+1, i+2)
									  + VAL2D(cur, j+1, i);	  
			}
		}
	
	#ifdef VERBOSE
	printf("Laplacien calcule\n");
	#endif

	ecrire(laplacien, (width-2), (height-2), file, "/laplacien");
	free(laplacien);
}

void derivee(double* cur, double* old, int width, int height, const char* file)
{
	int xx, yy;
	double *derivee = (double*)malloc(width*height*sizeof(double));
	
	for(xx=0; xx<width; xx++){
		for(yy=0; yy<height; yy++){
			VAL2D(derivee, xx, yy) = (VAL2D(old, xx, yy) - VAL2D(cur, xx, yy));
		}
	}
	
	#ifdef VERBOSE
	printf("Derivee calculee\n");
	#endif

	ecrire(derivee, width, height, file, "/Derivees");
	free(derivee);
}

// Ecriture dans le fichier HDF5

void ecriture(MPI_Comm cart_com, double* cur, int width, int height, const char* file, const char* set)
{
	int p, rang;
	MPI_Comm_rank(MPI_COMM_WORLD, &rang);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	int px, py;
	px=py=sqrt(p);
	
	int car_coord[2]; MPI_Cart_coords(cart_com, rang, 2, car_coord);
	
	hsize_t dim[2];
	hsize_t dim_g[2];
	dim[0]=height;
	dim_g[0]=(height-2)*px;
	dim[1]=width;
	dim_g[1]=(width-2)*py;
	
	hid_t file_id=-1;
	
	//Ouverture du fichier
	hid_t fapl_id=H5Pcreate(H5P_FILE_ACCESS);
	
	MPI_Info infos;
	MPI_Comm_get_info(cart_com, &infos);
	
	H5Pset_fapl_mpio(fapl_id, cart_com, infos);
	file_id=H5Fopen(file, H5F_ACC_RDWR, fapl_id);
	
	if(file_id < 0)
	{
		fprintf(stderr, "Impossible d'ouvrir le fichier\n");
		exit(EXIT_FAILURE);
	}
	else
		fprintf(stderr, "Ouverture reussie\n");
	
	// Creation du dataspace
	hid_t dataspace=-1;
	dataspace=H5Screate_simple(2, dim_g, dim_g);
	
	// Creation du dataset
	hid_t set_id=-1;
	set_id=H5Dcreate2(file_id, set, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	
	// Hyperslab
	hsize_t start[2];
	start[0]=(height-2)*car_coord[0];
	start[1]=(width-2)*car_coord[1];
	
	hsize_t stride[2];
	stride[0]=1;
	stride[1]=1;
	
	hsize_t count[2];
	count[0]=1;
	count[1]=1;
	
	hsize_t size[2];
	size[0]=height-2;
	size[1]=width-2;
	
	hid_t dataspace_local=-1;
	dataspace_local=H5Screate_simple(2, dim_g, dim_g);
	H5Sselect_hyperslab(dataspace_local, H5S_SELECT_SET, start, stride, size, NULL);
	
	//Ecriture
	H5Dwrite(set_id, H5T_NATIVE_DOUBLE, dataspace_local, dataspace, H5P_DEFAULT, cur);
	
	//Fermetures
	H5Pclose(fapl_id);
	H5Dclose(set_id);
	H5Sclose(dataspace_local);
	H5Sclose(dataspace);
	
	H5LTmake_dataset_double(file_id, set, 2, dim, cur);
	H5Fclose(file_id);
}

void init(double* dat, int width, int height, int px, int py)
{
	int xx, yy;
    py = py; // prevent unused warning
    for (yy=0; yy<height; ++yy) {
        for(xx=0; xx<width; ++xx) {
            VAL2D(dat,xx,yy) = 0;
        }
    }
    if ( px == 0 ) {
        for (yy=0; yy<height; ++yy) {
            VAL2D(dat,0,yy) = 1000000;
        }
    }
}

void iter(double* cur, double* next, int width, int height)
{
    int xx, yy;
    
    for(xx=0; xx<width; ++xx) {
        VAL2D(next,xx,0) = VAL2D(cur,xx,0);
    }
    for (yy=1; yy<height-1; ++yy) {
        VAL2D(next,0,yy) = VAL2D(cur,0,yy);
        for(xx=1; xx<width-1; ++xx) {
            VAL2D(next,xx,yy) =
                      (VAL2D(cur,xx,yy)   *.5)
                    + (VAL2D(cur,xx-1,yy) *.125)
                    + (VAL2D(cur,xx+1,yy) *.125)
                    + (VAL2D(cur,xx,yy-1) *.125)
                    + (VAL2D(cur,xx,yy+1) *.125);
        }
        VAL2D(next,width-1,yy) = VAL2D(cur,width-1,yy);
    }
    for(xx=0; xx<width; ++xx) {
        VAL2D(next,xx,height-1) = VAL2D(cur,xx,height-1);
    }

}

void exchange(MPI_Comm cart_com, double *cur, int width, int height)
{
    MPI_Status status;
    int rank_source, rank_dest;
    static MPI_Datatype column, row;
    static int initialized = 0;
   
    if ( !initialized ) {
        MPI_Type_vector(height-2, 1, width, MPI_DOUBLE, &column);
        MPI_Type_commit(&column);
        MPI_Type_contiguous(width-2, MPI_DOUBLE, &row);
        MPI_Type_commit(&row);
        initialized = 1;
    }
   
   
    /* send to the right */
    MPI_Cart_shift(cart_com, 0, 1, &rank_source, &rank_dest);
    MPI_Sendrecv(&VAL2D(cur, width-2, 1), 1, column, rank_dest,   100, /* send column before ghost */
                 &VAL2D(cur, 0,       1), 1, column, rank_source, 100, /* receive 1st column (ghost) */
                 cart_com, &status);

    /* send to the left */
    MPI_Cart_shift(cart_com, 0, -1, &rank_source, &rank_dest);
    MPI_Sendrecv(&VAL2D(cur, 1,       1), 1, column, rank_dest,   100, /* send column after ghost */
                 &VAL2D(cur, width-1, 1), 1, column, rank_source, 100, /* receive last column (ghost) */
                 cart_com, &status);

    /* send down */
    MPI_Cart_shift(cart_com, 1, 1, &rank_source, &rank_dest);
    MPI_Sendrecv(&VAL2D(cur, 1, height-2), 1, row, rank_dest,   100, /* send row before ghost */
                 &VAL2D(cur, 1, 0       ), 1, row, rank_source, 100, /* receive 1st row (ghost) */
                 cart_com, &status);

    /* send up */
    MPI_Cart_shift(cart_com, 1, -1, &rank_source, &rank_dest);
    MPI_Sendrecv(&VAL2D(cur, 1, 1       ), 1, row, rank_dest,   100, /* send column after ghost */
                 &VAL2D(cur, 1, height-1), 1, row, rank_source, 100, /* receive last column (ghost) */
                 cart_com, &status);
}

int main(int argc, char *argv[])
{
    //hid_t       file_id, data_set;   /* file identifier */
    //herr_t      status;
    MPI_Init(&argc, &argv);
   
    if ( argc != 4 ) {
        printf("Usage: %s <Nb_iter> <height> <width>\n", argv[0]);
        exit(1);
    }
   
    // Clean the file
	hid_t file_id=H5Fcreate("dataset.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	// hett_t status = H5LTmake_dataset(file_id, "data", RANK, dims, H5T_NATIVE_DOUBLE, cur);
	
	H5Fclose(file_id);
   
    int size; MPI_Comm_size(MPI_COMM_WORLD, &size);
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   
    int pwidth = sqrt(size);
    if ( ! pwidth ) pwidth=1;
    int pheight = size/pwidth;
    assert(pwidth * pheight == size);
   
    int nb_iter = atoi(argv[1]);
    int height = atoi(argv[2]);
    int width = atoi(argv[3]);
   
    //hsize_t dims[RANK]={height,width};
    // get local & add ghosts to sizes
    assert(width %pwidth ==0); width  = width /pwidth  + 2;
    assert(height%pheight==0); height = height/pheight + 2;
    assert(pwidth*pheight == size);
    
    int cart_dims[2] = { pwidth, pheight };
    int cart_period[2] = { 0, 0 };
    MPI_Comm cart_com; MPI_Cart_create(MPI_COMM_WORLD, 2, cart_dims, cart_period, 1, &cart_com);
    int car_coord[2]; MPI_Cart_coords(cart_com, rank, 2, car_coord);

    double *cur = malloc(sizeof(double)*width*height);
    double *next = malloc(sizeof(double)*width*height);
   
    init(cur, width, height, car_coord[0], car_coord[1]);

    int ii;
    for(ii=0; ii<nb_iter; ++ii) {
        iter(cur, next, width, height);
        exchange(cart_com, next, width, height);
        double *tmp = cur; cur = next; next = tmp;
    }
    
    ecriture(cart_com, cur, width, height, "dataset.h5", "/data");
	
	height = atoi(argv[2]);
    width  = atoi(argv[2]);
	lire(cur, width, height, "dataset.h5", "/data");
	moyennes(cur, width, height, size, rank,"dataset.h5");
	laplacien(cur, width, height, "dataset.h5");
	derivee(cur, next, width, height, "dataset.h5");

	free(cur);
	free(next);
    
    MPI_Finalize();
    return 0;
}
