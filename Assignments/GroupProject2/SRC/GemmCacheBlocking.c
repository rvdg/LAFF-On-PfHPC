#define A( i,j ) *( ap + (j)*lda + i )   // map A( i,j ) to array ap in column-major order
#define B( i,j ) *( bp + (j)*ldb + i )   // map B( i,j ) to array bp in column-major order
#define C( i,j ) *( cp + (j)*ldc + i )   // map C( i,j ) to array bp in column-major order

#define M_UNROLL 4
#define N_UNROLL 4

#define M_BLOCK 72
#define N_BLOCK 4080
#define K_BLOCK 256

#define min( i, j ) ( (i)<(j) ? (i): (j) )


#include<immintrin.h>


void GemmCacheBlockingMicroKernel(int k,
	      double *ap, int lda,
      	      double *bp, int ldb,
      	      double *cp, int ldc )
{

  const double* A_ptr = ap;
  const double* B_ptr = bp;
  double* C_ptr = cp;

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
  
  C_0123_0 = _mm256_add_pd(C_0123_0, _mm256_loadu_pd(C_ptr));
  C_0123_1 = _mm256_add_pd(C_0123_1, _mm256_loadu_pd(C_ptr + 1*ldc));
  C_0123_2 = _mm256_add_pd(C_0123_2, _mm256_loadu_pd(C_ptr + 2*ldc));
  C_0123_3 = _mm256_add_pd(C_0123_3, _mm256_loadu_pd(C_ptr + 3*ldc));

  _mm256_storeu_pd(C_ptr        , C_0123_0);
  _mm256_storeu_pd(C_ptr + 1*ldc, C_0123_1);
  _mm256_storeu_pd(C_ptr + 2*ldc, C_0123_2);
  _mm256_storeu_pd(C_ptr + 3*ldc, C_0123_3);

}

void GemmCacheBlockingInnerKernel( int m, int n, int k,
	      double *ap, int lda,
      	      double *bp, int ldb,
      	      double *cp, int ldc )
{
  int i, j;

  double *a_part, *b_part, *c_part;

  for ( j = 0;j < n;j += N_UNROLL)
  {
    b_part = bp + (j)*ldb;
    for ( i = 0;i < m;i += M_UNROLL)
    {
      a_part = ap + i;
      c_part = cp + j*ldc + i;
      GemmCacheBlockingMicroKernel( k, a_part, lda, b_part, ldb, c_part, ldc );
    }
  }

  return;
}
  

void GemmCacheBlocking( int m, int n, int k,
	      double *ap, int lda,
      	      double *bp, int ldb,
      	      double *cp, int ldc )
{
  int i, j, p;

  double *b_part, *a_part, *c_part;

  for (j = 0;j < n;j += N_BLOCK)
  {
    int n_part = min(N_BLOCK, n-j);

    for (p = 0;p < k;p += K_BLOCK)
    {
      int k_part = min(K_BLOCK, k-p);

      b_part = &B(p, j);

      for (i = 0;i < m;i += M_BLOCK)
      {
        int m_part = min(M_BLOCK, m-i);

        a_part = &A(i, p);
        c_part = &C(i, j);

        GemmCacheBlockingInnerKernel(m_part, n_part, k_part, a_part, lda, b_part, ldb, c_part, ldc);
      }
    }
  }

  return;
}
  
