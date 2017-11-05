#include <stdio.h>
#include <stdlib.h>

#define alpha( i,j ) A[ (j)*ldA + (i) ]   // map alpha( i,j ) to array A
#define beta( i,j )  B[ (j)*ldB + (i) ]   // map beta( i,j ) to array B
#define gamma( i,j ) C[ (j)*ldC + (i) ]   // map gamma( i,j ) to array C

#define min( x, y ) ( ( x ) < ( y ) ? x : y )

#define NC 1024
#define KC 128
#define MC 128
#define MR 4
#define NR 4

void GemmLoopFive( int, int, int, double *, int, double *, int, double *, int );
void GemmLoopFour( int, int, int, double *, int, double *, int, double *, int );
void GemmLoopThree( int, int, int, double *, int, double *, double *, int );
void GemmLoopTwo( int, int, int, double *, double *, double *, int );
void GemmLoopOne( int, int, int, double *, double *, double *, int );
void MicroKernel_MRxNR( int, int, int, double *, double *, double *, int );
void PackBlockA_MCxKC( int, int, double *, int, double * );
void PackPanelB_KCxNC( int, int, double *, int, double * );
  
void GemmFiveLoops( int m, int n, int k, double *A, int ldA,
		   double *B, int ldB, double *C, int ldC )
{
  if ( m % MR != 0 || MC % MR != 0 ){
    printf( "m and MC must be multiples of MR\n" );
    exit( 0 );
  }
  if ( n % NR != 0 || NC % MR != 0 ){
    printf( "n and NC must be multiples of NR\n" );
    exit( 0 );
  }

  GemmLoopFive( m, n, k, A, ldA, B, ldB, C, ldC );
}

void GemmLoopFive( int m, int n, int k, double *A, int ldA,
		   double *B, int ldB, double *C, int ldC )
{
  for ( int j=0; j<n; j+=NC ) {
    int jb = min( NC, n-j );    /* Last loop may not involve a full block */

    GemmLoopFour( m, jb, k, A, ldA, &beta( 0,j ), ldB, &gamma( 0,j ), ldC );
  }
}

void GemmLoopFour( int m, int n, int k, double *A, int ldA,
		   double *B, int ldB, double *C, int ldC )
{
  double *PackedPanelB = ( double * ) malloc( KC * NC * sizeof( double ) );
  
  for ( int p=0; p<k; p+=KC ) {
    int pb = min( KC, k-p );    /* Last loop may not involve a full block */

    PackPanelB_KCxNC( pb, n, &beta( p, 0 ), ldB, PackedPanelB );
    GemmLoopThree( m, n, pb, &alpha( 0, p ), ldA, PackedPanelB, C, ldC );
  }
}

void GemmLoopThree( int m, int n, int k, double *A, int ldA,
		   double *PackedPanelB, double *C, int ldC )
{
  double *PackedBlockA = ( double * ) malloc( MC * KC * sizeof( double ) );
  
  for ( int i=0; i<m; i+=MC ) {
    int ib = min( MC, m-i );    /* Last loop may not involve a full block */
    
    PackBlockA_MCxKC( ib, k, &alpha( i, 0 ), ldA, PackedBlockA );

    GemmLoopTwo( ib, n, k, PackedBlockA, PackedPanelB, &gamma( i,0 ), ldC );
  }
}

void GemmLoopTwo( int m, int n, int k, double *PackedBlockA, 
		  double *PackedPanelB, double *C, int ldC )
{
  for ( int j=0; j<n; j+=NR ) {
    int jb = min( NR, n-j );

    GemmLoopOne( m, jb, k, PackedBlockA, &PackedPanelB[ j*k ], &gamma( 0,j ), ldC );
  }
}

void GemmLoopOne( int m, int n, int k, double *PackedBlockA, 
		  double *MicroPanelB, double *C, int ldC )
{
  for ( int i=0; i<m; i+=MR ) {
    int ib = min( MR, m-i );
    
    MicroKernel_MRxNR( ib, n, k, &PackedBlockA[ i*k ], MicroPanelB, &gamma( i,0 ), ldC );
  }
}

