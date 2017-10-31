#define alpha( i,j ) A[ (j)*ldA + i ]   // map alpha( i,j ) to array A 
#define beta( i,j )  B[ (j)*ldB + i ]   // map beta( i,j )  to array B
#define gamma( i,j ) C[ (j)*ldC + i ]   // map gamma( i,j ) to array C

void GemvIJ_Dots( int, int, double *, int, double *, int, double *, int );

void GemmJIP_Gemv( int m, int n, int k,
		   double *A, int ldA,
		   double *B, int ldB,
		   double *C, int ldC )
{
  int j;

  for ( j=0; j<n; j++ )
    GemvIJ_Dots( m, k, A, ldA, &beta( 0,j ), 1, &gamma( 0,j ), 1 ); 

  return;
}
  
