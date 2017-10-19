#define A( i,j ) *( ap + (j)*lda + i )   // map A( i,j ) to array ap in column-major order
#define B( i,j ) *( bp + (j)*ldb + i )   // map B( i,j ) to array bp in column-major order
#define C( i,j ) *( cp + (j)*ldc + i )   // map C( i,j ) to array bp in column-major order

#define M_UNROLL 4
#define N_UNROLL 4


void GemmUnrollKernel(int k,
	      double *ap, int lda,
      	      double *bp, int ldb,
      	      double *cp, int ldc )
{
    for (int p = 0;p < k;p++)
    {
        C(0,0) += A(0,p) * B(p,0); 
	C(0,1) += A(0,p) * B(p,1);
        C(0,2) += A(0,p) * B(p,2); 
	C(0,3) += A(0,p) * B(p,3);

        C(1,0) += A(1,p) * B(p,0); 
	C(1,1) += A(1,p) * B(p,1);
        C(1,2) += A(1,p) * B(p,2); 
	C(1,3) += A(1,p) * B(p,3);

        C(2,0) += A(2,p) * B(p,0); 
	C(2,1) += A(2,p) * B(p,1);
        C(2,2) += A(2,p) * B(p,2); 
	C(2,3) += A(2,p) * B(p,3);     

        C(3,0) += A(3,p) * B(p,0); 
	C(3,1) += A(3,p) * B(p,1);
        C(3,2) += A(3,p) * B(p,2); 
	C(3,3) += A(3,p) * B(p,3);
        
    }

}

void GemmUnroll( int m, int n, int k,
	      double *ap, int lda,
      	      double *bp, int ldb,
      	      double *cp, int ldc )
{
  int i, j;

  for ( j = 0;j < n;j += N_UNROLL)
  {
    for ( i = 0;i < m;i += M_UNROLL)
    {
      GemmUnrollKernel( k, &A( i,0 ), lda, &B( 0, j ), ldb, &C(i,j ), ldc );
    }
  }

  return;
}
  
