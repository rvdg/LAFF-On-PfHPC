#define alpha( i,j ) A[ (j)*ldA + (i) ]   // map alpha( i,j ) to array A
#define beta( i,j )  B[ (j)*ldB + (i) ]   // map beta( i,j ) to array B
#define gamma( i,j ) C[ (j)*ldC + (i) ]   // map gamma( i,j ) to array C

#define NC 1024
#define KC 128
#define MC 128
#define MR 4
#define NR 4

void GemmLoopFive( int, int, int, double *, int, double *, int, double *, int );
void GemmLoopFour( int, int, int, double *, int, double *, int, double *, int );
void GemmLoopThree( int, int, int, double *, int, double *, int, double *, int );
void GemmLoopTwo( int, int, int, double *, int, double *, int, double *, int );
void GemmLoopOne( int, int, int, double *, int, double *, int, double *, int );
void GemmIntrinsicsKernel_mrxnr( int, double *, int, double *, int,
				      double *, int );

void GemmLoopFive( int m, int n, int k, double *A, int ldA,
		   double *B, int ldB, double *C, int ldC )
{
  int nb;

  if ( m % MR ~= 0 || MC % MR ~= 0 ){
    printf( "m and MC must be multiples of MR\n" );
    exit( 0 );
  }
  if ( n % NR ~= 0 || NC % MR ~= 0 ){
    printf( "n and NC must be multiples of NR\n" );
    exit( 0 );
  }
    
  for ( int j=0; j<n; j+=NC ) {
    nb = min( NC, n-j );    /* Last loop may not involve a full block */
    GemmLoopFour( m, jb, k, A, ldA, &beta( 0,j ), ldB, &gamma( 0,j ), ldC );
  }
}

void GemmLoopFour( int m, int n, int k, double *A, int ldA,
		   double *B, int ldB, double *C, int ldC )
{
  int kb;
  
  for ( int p=0; p<k; p+=KC ) {
    kb = min( KC, k-p );    /* Last loop may not involve a full block */
    GemmLoopThree( m, n, kb, &alpha( 0, p ), ldA, &beta( p,0 ), ldB, C, ldC );
  }
}

void GemmLoopThree( int m, int n, int k, double *A, int ldA,
		   double *B, int ldB, double *C, int ldC )
{
  int mb;
  
  for ( int i=0; i<m; i+=MC ) {
    mb = min( MC, m-i );    /* Last loop may not involve a full block */
    GemmLoopTwo( mb, n, k, &alpha( i,0 ), ldA, B, ldB, &gamma( i,0 ), ldC );
  }
}

void GemmLoopTwo( int m, int n, int k, double *A, int ldA,
		  double *B, int ldB, double *C, int ldC )
{
  int nb;
  
  for ( int j=0; j<n; j+=NR ) {
    nb = min( NR, n-j );    
    GemmLoopOne( m, nb, k, A, ldA, &beta( 0,j ), ldB, &gamma( 0,j ), ldC );
  }
}

void GemmLoopOne( int m, int n, int k, double *A, int ldA,
		  double *B, int ldB, double *C, int ldC )
{
  int mb;
  
  for ( int i=0; i<m; i+=MR ) {
    mb = min( MR, m-i );
    GemmIntrinsicsKernel_mrxnr( k, &alpha( i,0 ), ldA, B, ldB, &gamma( i,0 ), ldC );
  }
}

#include<immintrin.h>

void GemmIntrinsicsKernel_mrxnr( int k, double *A, int ldA,
				 double *B, int ldB, double *C, int ldC )
{
  __m256d gamma_0123_0 = _mm256_loadu_pd( &gamma( 0,0 ) );
  __m256d gamma_0123_1 = _mm256_loadu_pd( &gamma( 0,1 ) );
  __m256d gamma_0123_2 = _mm256_loadu_pd( &gamma( 0,2 ) );
  __m256d gamma_0123_3 = _mm256_loadu_pd( &gamma( 0,3 ) );
   	
  for ( int p=0; p<k; p++ ){
    /* load alpha( 0:3, p ) */
    __m256d alpha_0123_p = _mm256_loadu_pd( &alpha( 0,p ) );
    /* load beta( p, 0 ); update gamma( 0:3, 0 ) */
    __m256d beta_p_j = _mm256_broadcast_sd( &beta( p, 0) );
    gamma_0123_0 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_0 );
    /* load beta( p, 1 ); update gamma( 0:3, 1 ) */
    __m256d beta_p_j = _mm256_broadcast_sd( &beta( p, 1) );
    gamma_0123_1 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_1 );
    /* load beta( p, 2 ); update gamma( 0:3, 2 ) */
    __m256d beta_p_j = _mm256_broadcast_sd( &beta( p, 2) );
    gamma_0123_2 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_2 );
    /* load beta( p, 3 ); update gamma( 0:3, 3 ) */
    __m256d beta_p_j = _mm256_broadcast_sd( &beta( p, 3) );
    gamma_0123_3 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_3 ); 
  }

  /* Store the updated results */
  _mm256_storeu_pd( &gamma(0,0), gamma_0123_0 );
  _mm256_storeu_pd( &gamma(0,1), gamma_0123_1 );
  _mm256_storeu_pd( &gamma(0,2), gamma_0123_2 );
  _mm256_storeu_pd( &gamma(0,3), gamma_0123_3 );

}
