#define A( i,j ) *( ap + (j)*lda + i )   // map A( i,j ) to array ap in column-major order
#define B( i,j ) *( bp + (j)*ldb + i )   // map B( i,j ) to array bp in column-major order
#define C( i,j ) *( cp + (j)*ldc + i )   // map C( i,j ) to array bp in column-major order


void GemmIJP( int m, int n, int k,
	      double *ap, int lda,
      	      double *bp, int ldb,
      	      double *cp, int ldc )
{
  int i, j, p;

  for ( i=0; i<m; i++ )
    for ( j=0; j<n; j++ )
      for ( p=0; p<k; p++ )
	C( i,j ) = A( i,p ) * B( p, j ) + C( i,j );

  return;
}
  
