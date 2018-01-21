#include <stdio.h>
#include <stdlib.h>

#define alpha( i,j ) A[ (j)*ldA + i ]   // map alpha( i,j ) to array A 
#define beta( i,j )  B[ (j)*ldB + i ]   // map beta( i,j )  to array B
#define gamma( i,j ) C[ (j)*ldC + i ]   // map gamma( i,j ) to array C

void Gemm_J_Gemv( int, int, int, double *, int, double *, int, double *, int );
void dgemv_( char *, int *, int *, double *, double *, int *,   double *, int *,
                                   double *, double *, int * );

void GemmWrapper( int m, int n, int k, double *A, int ldA,
                                       double *B, int ldB, double *C, int ldC )
{
  Gemm_J_Gemv( m, n, k, A, ldA, B, ldB, C, ldC );
}

void Gemm_J_Gemv( int m, int n, int k, double *A, int ldA, double *B, int ldB,
                                       double *C, int ldC )
{
  int i_one=1;
  double d_one=1.0;

  for ( int j=0; j<n; j++ )
    dgemv_( "N", &m, &k, &d_one, A, &ldA, &beta( 0,j ), &i_one,
                         &d_one, &gamma( 0, j ), &i_one );
}
  
