#define A( i,j ) *( ap + (j)*lda + i )   // map A( i,j ) to array ap in column-major order
#define B( i,j ) *( bp + (j)*ldb + i )   // map B( i,j ) to array bp in column-major order
#define C( i,j ) *( cp + (j)*ldc + i )   // map C( i,j ) to array bp in column-major order

void Sapdot( int, double *, int, double *, double * );
  
void GemmIJSapdot( int m, int n, int k,
	   	   double *ap, int lda,
		   double *bp, int ldb,
		   double *cp, int ldc )
{
  int i, j;

  for ( i=0; i<m; i++ )
    for ( j=0; j<n; j++ )
      Sapdot( k, &A( i,0 ), lda, &B( 0, j ), &C(i,j ) );

  return;
}

#define X( i )  x[ (i)*incx ]

void Sapdot( int n,
	     double *x, int incx,
	     double *y, 
	     double *gamma )
{
  int i;
  
  for ( i=0; i<n; i++ )
    *gamma += X( i ) * y[ i ];
}
    
