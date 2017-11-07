#include <stdio.h>
#include <stdlib.h>
#include<immintrin.h>

#define alpha( i,j ) A[ (j)*ldA + (i) ]   // map alpha( i,j ) to array A
#define beta( i,j )  B[ (j)*ldB + (i) ]   // map beta( i,j ) to array B
#define gamma( i,j ) C[ (j)*ldC + (i) ]   // map gamma( i,j ) to array C

#define MR 4
#define NR 4

#define MC 128
#define NC 128
#define KC 128

#define min( x, y )  ( (x) < (y) ? x : y )

void GemmIntrinsicsKernel_MRxNR( int, double *, int, double *, int,
		      double *, int );

void GemmJI_Kernel_MRxNR( int, int, int, double *, int, double *, int,
		      double *, int );

void GemmIJP_JI_Kernel_MRxNR( int, int, int, double *, int, double *, int,
		      double *, int );

void GemmWrapper( int m, int n, int k, double *A, int ldA,
		  double *B, int ldB, double *C, int ldC )
{
  if ( m % MR != 0 || MC % MR != 0 ){
    printf( "m and MC must be multiples of MR\n" );
    exit( 0 );
  }
  if ( n % NR != 0 || NC % NR != 0 ){
    printf( "n and NC must be multiples of NR\n" );
    exit( 0 );
  }

  GemmIJP_JI_Kernel_MRxNR( m, n, k, A, ldA, B, ldB, C, ldC );
}

void GemmIJP_JI_Kernel_MRxNR( int m, int n, int k, double *A, int ldA,
			  double *B, int ldB, double *C, int ldC )
{
  int ib, jb, pb;
  
  for ( int i=0; i<m; i+=MC ) {
    ib = min( MC, m-i );        /* Last block may not be a full block */
    for ( int j=0; j<n; j+=NC ) {
      jb = min( NC, n-j );        /* Last block may not be a full block */
      for ( int p=0; p<k; p+=KC ) {
	pb = min( KC, k-p );        /* Last block may not be a full block */

	GemmJI_Kernel_MRxNR
	  ( ib, jb, pb, &alpha( i,p ), ldA, &beta( p,j ), ldB, &gamma( i,j ), ldC );
      }
    }
  }
}

void GemmJI_Kernel_MRxNR( int m, int n, int k, double *A, int ldA,
			  double *B, int ldB, double *C, int ldC )
{
  for ( int j=0; j<n; j+=NR ) /* n is assumed to be a multiple of NR */
    for ( int i=0; i<m; i+=MR ) /* m is assumed to be a multiple of MR */
      GemmIntrinsicsKernel_MRxNR
        ( k, &alpha( i,0 ), ldA, &beta( 0,j ), ldB, &gamma( i,j ), ldC );
}


void GemmIntrinsicsKernel_MRxNR( int k, double *A, int ldA,
				 double *B, int ldB, double *C, int ldC )
{
  __m256d gamma_0123_0 = _mm256_loadu_pd( &gamma( 0,0 ) );
  __m256d gamma_0123_1 = _mm256_loadu_pd( &gamma( 0,1 ) );
  __m256d gamma_0123_2 = _mm256_loadu_pd( &gamma( 0,2 ) );
  __m256d gamma_0123_3 = _mm256_loadu_pd( &gamma( 0,3 ) );

  __m256d beta_p_j;
   	
  for ( int p=0; p<k; p++ ){
    /* load alpha( 0:3, p ) */
    __m256d alpha_0123_p = _mm256_loadu_pd( &alpha( 0,p ) );
    /* load beta( p, 0 ); update gamma( 0:3, 0 ) */
    beta_p_j = _mm256_broadcast_sd( &beta( p,0) );
    gamma_0123_0 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_0 );
    /* load beta( p, 1 ); update gamma( 0:3, 1 ) */
    beta_p_j = _mm256_broadcast_sd( &beta( p,1) );
    gamma_0123_1 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_1 );
    /* load beta( p, 2 ); update gamma( 0:3, 2 ) */
    beta_p_j = _mm256_broadcast_sd( &beta( p,2) );
    gamma_0123_2 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_2 );
    /* load beta( p, 3 ); update gamma( 0:3, 3 ) */
    beta_p_j = _mm256_broadcast_sd( &beta( p,3) );
    gamma_0123_3 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_3 ); 
  }

  /* Store the updated results */
  _mm256_storeu_pd( &gamma(0,0), gamma_0123_0 );
  _mm256_storeu_pd( &gamma(0,1), gamma_0123_1 );
  _mm256_storeu_pd( &gamma(0,2), gamma_0123_2 );
  _mm256_storeu_pd( &gamma(0,3), gamma_0123_3 );
}
