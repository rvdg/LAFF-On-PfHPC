#define alpha( i,j ) A[ (j)*ldA + i ]   // map alpha( i,j ) to array A 
#define chi( i )   x[ (i) * incx ]      // map chi( i ) to array x

void Axpy( int, double, double *, int, double *, int );

void GemvJI_Axpy( int m, int n, double *A, int ldA, double *x, int incx,
                                double *y, int incy )
{
  int j;

  for ( j=0; j<n; j++ )
    Axpy( m, chi( j ), &alpha( 0,j ), 1, y, incy );

  return;
}
  
