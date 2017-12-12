#define A( i,j ) a[ (j)*lda + i ]   // map A( i,j ) to array a 
#define B( i,j ) b[ (j)*ldb + i ]   // map B( i,j ) to array b
#define C( i,j ) c[ (j)*ldc + i ]   // map C( i,j ) to array c

void Dots( int, double *, int, double *, int, double * );

void GemmIJP_Dots( int m, int n, int k,
	      double *a, int lda,
	      double *b, int ldb,
	      double *c, int ldc )
{
  int i, j;

  for ( i=0; i<m; i++ )
    for ( j=0; j<n; j++ )
      Dots( k, &A( i,0 ), lda, &B( 0,j ), 1, &C( i,j ) );

  return;
}
  
