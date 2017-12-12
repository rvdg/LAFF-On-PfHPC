#include <stdio.h>
#include <stdlib.h>

#define alpha( i,j ) A[ (j)*ldA + i ]   // map alpha( i,j ) to array A 
#define beta( i,j )  B[ (j)*ldB + i ]   // map beta( i,j )  to array B
#define gamma( i,j ) C[ (j)*ldC + i ]   // map gamma( i,j ) to array C

void Gemm_IPJ( int, int, int, double *, int, double *, int, double *, int );

void GemmWrapper( int m, int n, int k, double *A, int ldA,
		  double *B, int ldB, double *C, int ldC )
{
  Gemm_IPJ( m, n, k, A, ldA, B, ldB, C, ldC );
}


void Gemm_IPJ( int m, int n, int k,
               double *A, int ldA,
               double *B, int ldB,
               double *C, int ldC )
{
  int i, j, p;

  for ( i=0; i<m; i++ )
    for ( p=0; p<k; p++ )
      for ( j=0; j<n; j++ )
        gamma( i,j ) += alpha( i,p ) * beta( p,j );

  return;
}
  
