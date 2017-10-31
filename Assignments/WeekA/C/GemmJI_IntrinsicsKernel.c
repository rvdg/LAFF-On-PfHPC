#define alpha( i,j ) A[ (j)*ldA + (i) ]   // map alpha( i,j ) to array A
#define beta( i,j )  B[ (j)*ldB + (i) ]   // map beta( i,j ) to array B
#define gamma( i,j ) C[ (j)*ldC + (i) ]   // map gamma( i,j ) to array C

#define MR 4
#define NR 4

void GemmIntrinsicsKernel_m4xn4( int  , double *, int ,
                                        double *, int ,
                                        double *, int );

void GemmJI_IntrinsicsKernel( int m, int n, int k,
                              double *A, int ldA,
                              double *B, int ldB,
                              double *C, int ldC )
{
  int i, j;

  for ( int j=0; j<n; j+=NR ) /* n is assumed to be a multiple of NR */
    for ( int i=0; i<m; i+=MR ) /* m is assumed to be a multiple of MR */
      GemmIntrinsicsKernel_m4xn4
        ( k, &alpha( i,0 ), ldA, &beta( 0,j ), ldB, &gamma( i,j ), ldC );
  
  return;
}
