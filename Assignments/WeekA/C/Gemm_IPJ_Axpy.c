#define alpha( i,j ) A[ (j)*ldA + i ]   // map alpha( i,j ) to array A 
#define beta( i,j )  B[ (j)*ldB + i ]   // map beta( i,j )  to array B
#define gamma( i,j ) C[ (j)*ldC + i ]   // map gamma( i,j ) to array C

void Axpy( int n, double alpha, double *x, int incx, double *y, int incy );

void GemmIPJ_Axpy( int m, int n, int k,
		   double *A, int ldA,
		   double *B, int ldB,
		   double *C, int ldC )
{
  int i, p;

  for ( i=0; i<m; i++ )
    for ( p=0; p<k; p++ )
       Axpy( n,    ,    ,    ,    ,    );

  return;
}
  
