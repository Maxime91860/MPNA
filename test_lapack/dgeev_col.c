#include <stdlib.h>
#include <stdio.h>
#include <lapacke.h>

/* Parameters */
#define N 3
#define LDA N
#define LDVL N
#define LDVR N

/* Auxiliary routine: printing eigenvalues */
void print_eigenvalues( char* desc, int n, double* wr, double* wi ) {
        int j;
        printf( "\n %s\n", desc );
   for( j = 0; j < n; j++ ) {
      if( wi[j] == (double)0.0 ) {
         printf( " %6.2f", wr[j] );
      } else {
         printf( " (%6.2f,%6.2f)", wr[j], wi[j] );
      }
   }
   printf( "\n" );
}

/* Auxiliary routine: printing eigenvectors */
void print_eigenvectors( char* desc, int n, double* wi, double* v, int ldv ) {
        int i, j;
        printf( "\n %s\n", desc );
   for( i = 0; i < n; i++ ) {
      j = 0;
      while( j < n ) {
         if( wi[j] == (double)0.0 ) {
            printf( " %6.2f", v[i+j*ldv] );
            j++;
         } else {
            printf( " (%6.2f,%6.2f)", v[i+j*ldv], v[i+(j+1)*ldv] );
            printf( " (%6.2f,%6.2f)", v[i+j*ldv], -v[i+(j+1)*ldv] );
            j += 2;
         }
      }
      printf( "\n" );
   }
}

int main(int argc, char const *argv[])
{
    int n = N, lda = LDA, ldvl = LDVL, ldvr = LDVR, info;
    /* Local arrays */
    double wr[N], wi[N], vl[LDVL*N], vr[LDVR*N];
    double a[LDA*N] = {
       3.00,  -1.00,  -3.00, 
       -1.00,  3.00,  -1.00, 
       0.00, 0.00, 6.00
    };	
    /* Executable statements */
    printf( "LAPACKE_dgeev (row-major, high-level) Example Program Results\n" );
    /* Solve eigenproblem */
    info = LAPACKE_dgeev( LAPACK_ROW_MAJOR, 'V', 'V', n, a, lda, wr, wi,
                    vl, ldvl, vr, ldvr );
    /* Check for convergence */
    if( info > 0 ) {
            printf( "The algorithm failed to compute eigenvalues.\n" );
            exit( 1 );
    }
        /* Print eigenvalues */
        print_eigenvalues( "Eigenvalues", n, wr, wi );
        /* Print left eigenvectors */
        print_eigenvectors( "Left eigenvectors", n, wi, vl, ldvl );
        /* Print right eigenvectors */
        print_eigenvectors( "Right eigenvectors", n, wi, vr, ldvr );
        exit( 0 );
} /* End of LAPACKE_dgeev Example */

