#define alpha( i,j ) A[ (j)*ldA + i ]   // map alpha( i,j ) to array A 
#define beta( i,j )  B[ (j)*ldB + i ]   // map beta( i,j )  to array B
#define gamma( i,j ) C[ (j)*ldC + i ]   // map gamma( i,j ) to array C

void GerJI_Axpy( int, int, double *, int, double *, int, double *, int );

void GemmPJI_Ger( int m, int n, int k,
		  double *A, int ldA,
		  double *B, int ldB,
		  double *C, int ldC )
{
  int p;

  for ( p=0; p<k; p++ )
    GerJI_Axpy(m, n, &alpha(0, p), 1, &beta(p, 0), ldB, C, ldC); 

  return;
}
  
