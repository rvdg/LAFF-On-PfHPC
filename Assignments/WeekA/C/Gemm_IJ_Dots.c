#include <stdio.h>
#include <stdlib.h>

#define alpha( i,j ) A[ (j)*ldA + i ]   // map alpha( i,j ) to array A 
#define beta( i,j )  B[ (j)*ldB + i ]   // map beta( i,j )  to array B
#define gamma( i,j ) C[ (j)*ldC + i ]   // map gamma( i,j ) to array C

void Gemm_IJ_Dots( int, int, int, double *, int, double *, int, double *, int );
void Dots( int, double *, int, double *, int, double * );
double ddot_( int *, double *, int *, double *, int * );

void GemmWrapper( int m, int n, int k, double *A, int ldA,
		  double *B, int ldB, double *C, int ldC )
{
  Gemm_IJ_Dots( m, n, k, A, ldA, B, ldB, C, ldC );
}


void Gemm_IJ_Dots( int m, int n, int k,
		   double *A, int ldA,
		   double *B, int ldB,
		   double *C, int ldC )
{
  int i, j, i_one=1;

  for ( i=0; i<m; i++ )
    for ( j=0; j<n; j++ )
      //      Dots( k, &alpha( i,0 ), ldA, &beta( 0,j ), 1, &gamma( i,j ) );
      gamma( i,j ) += ddot_( &k, &alpha( i,0 ), &ldA, &beta( 0,j ), &i_one );
}

#define chi( i ) x[ (i)*incx ]   // map chi( i ) to array x 
#define psi( i ) x[ (i)*incy ]   // map psi( i ) to array y

void Dots( int n, double *x, int incx, double *y, int incy, double *gamma )
{
  int i;

  for ( i=0; i<n; i++ )
    *gamma += chi( i ) * psi( i );
  
  return;
}

  
