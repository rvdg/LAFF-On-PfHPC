#define A( i,j ) *( ap + (j)*lda + i )   // map A( i,j ) to array ap in column-major order
#define B( i,j ) *( bp + (j)*ldb + i )   // map B( i,j ) to array bp in column-major order
#define C( i,j ) *( cp + (j)*ldc + i )   // map C( i,j ) to array bp in column-major order

#define M_UNROLL 4
#define N_UNROLL 4

#include<immintrin.h>


void GemmIntrinsicsKernel(int k,
	      double *ap, int lda,
      	      double *bp, int ldb,
      	      double *cp, int ldc )
{

  const double* A_ptr = ap;
  const double* B_ptr = bp;

    __m256d C_0123_0 = _mm256_setzero_pd();
    __m256d C_0123_1 = _mm256_setzero_pd();
    __m256d C_0123_2 = _mm256_setzero_pd();
    __m256d C_0123_3 = _mm256_setzero_pd();
   	
  for (int p = 0;p < k;p++)
  {
    __m256d A_0123 = _mm256_loadu_pd(A_ptr + 0);

    __m256d B_0 = _mm256_broadcast_sd(B_ptr + 0*ldb);
    __m256d B_1 = _mm256_broadcast_sd(B_ptr + 1*ldb);
    C_0123_0 = _mm256_add_pd(C_0123_0, _mm256_mul_pd(A_0123, B_0));
    C_0123_1 = _mm256_add_pd(C_0123_1, _mm256_mul_pd(A_0123, B_1));

    __m256d B_2 = _mm256_broadcast_sd(B_ptr + 2*ldb);
    __m256d B_3 = _mm256_broadcast_sd(B_ptr + 3*ldb);
    C_0123_2 = _mm256_add_pd(C_0123_0, _mm256_mul_pd(A_0123, B_2));
    C_0123_3 = _mm256_add_pd(C_0123_1, _mm256_mul_pd(A_0123, B_3));

     A_ptr += lda;
     B_ptr ++;
  }
  
  C_0123_0 = _mm256_add_pd(C_0123_0, _mm256_loadu_pd(&C(0,0)));
  C_0123_1 = _mm256_add_pd(C_0123_1, _mm256_loadu_pd(&C(0,1)));
  C_0123_2 = _mm256_add_pd(C_0123_2, _mm256_loadu_pd(&C(0,2)));
  C_0123_3 = _mm256_add_pd(C_0123_3, _mm256_loadu_pd(&C(0,3)));

  _mm256_storeu_pd(&C(0,0), C_0123_0);
  _mm256_storeu_pd(&C(0,1), C_0123_1);
  _mm256_storeu_pd(&C(0,2), C_0123_2);
  _mm256_storeu_pd(&C(0,3), C_0123_3);

}

void GemmIntrinsics( int m, int n, int k,
	      double *ap, int lda,
      	      double *bp, int ldb,
      	      double *cp, int ldc )
{
  int i, j;

  for ( j = 0;j < n;j += N_UNROLL)
  {
    for ( i = 0;i < m;i += M_UNROLL)
    {
      GemmIntrinsicsKernel( k, &A( i,0 ), lda, &B( 0, j ), ldb, &C(i,j ), ldc );
    }
  }

  return;
}
  
