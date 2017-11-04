#define alpha( i,j ) A[ (j)*ldA + (i) ]   // map alpha( i,j ) to array A
#define beta( i,j )  B[ (j)*ldB + (i) ]   // map beta( i,j ) to array B
#define gamma( i,j ) C[ (j)*ldC + (i) ]   // map gamma( i,j ) to array C

#define MR 4
#define NR 4

#define MC 128
#define NC 128
#define KC 128

void GemmIntrinsicsKernel_MRxNR( int, double *, int, double *, int,
		      double *, int );

void GemmJI_Kernel_MRxNR( int, int, int, double *, int, double *, int,
		      double *, int );

void GemmIJP_JI_Kernel_MRxNR( int m, int n, int k, double *A, int ldA,
			  double *B, int ldB, double *C, int ldC )
{
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

