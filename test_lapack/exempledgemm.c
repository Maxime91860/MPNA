

#include  <stdlib.h>
#include  <stdio.h>
#define  dgemm  dgemm_
#define  dlarnv  dlarnv_
#define  MATSIZE  256
static  int  IONE = 1;

static  double  c__1 = 1.0;
// char  NoTranspose = ’N’;
char NoTranspose = 'N';
int  ISEED [4] = {0,0,0,1};

int  main (int argc , char **argv){
	double *A, *B, *C;	
	int N = MATSIZE;
	int NN = N*N;
	/*  Allocate  data */
	A     = (double  *) malloc(N*N*sizeof(double ));
	B     = (double  *) malloc(N*N*sizeof(double ));
	C     = (double  *) malloc(N*N*sizeof(double ));

	/* Fill  the  matrices  with  random  values  */
	dlarnv( &IONE , ISEED , &NN , A );
	dlarnv( &IONE , ISEED , &NN , B );
	dlarnv( &IONE , ISEED , &NN , C );

	/*  Perform  the  operation: C = 1.0*A*B + 1.0*C */
	dgemm(&NoTranspose, &NoTranspose, &N, &N, &N, &c__1, A, &N, B, &N, &c__1 , C, &N);
	
	free(A);
	free(B);
	free(C);
	return  EXIT_SUCCESS;
}