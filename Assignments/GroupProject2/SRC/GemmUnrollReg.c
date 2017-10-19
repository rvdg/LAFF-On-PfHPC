#define A( i,j ) *( ap + (j)*lda + i )   // map A( i,j ) to array ap in column-major order
#define B( i,j ) *( bp + (j)*ldb + i )   // map B( i,j ) to array bp in column-major order
#define C( i,j ) *( cp + (j)*ldc + i )   // map C( i,j ) to array bp in column-major order

#define M_UNROLL 4
#define N_UNROLL 4


void GemmUnrollRegKernel(int k,
	      double *ap, int lda,
      	      double *bp, int ldb,
      	      double *cp, int ldc )
{

  const double* A_ptr = ap;
  const double* B_ptr = bp;

  double C_tmp[M_UNROLL][N_UNROLL] = {};
   	
  for (int p = 0;p < k;p++)
  {
     double B_value = B_ptr[0*ldb];
     C_tmp[0][0] += A_ptr[0] * B_value;
     C_tmp[1][0] += A_ptr[1] * B_value;
     C_tmp[2][0] += A_ptr[2] * B_value;
     C_tmp[3][0] += A_ptr[3] * B_value;
    
     B_value = B_ptr[1*ldb];
     C_tmp[0][1] += A_ptr[0] * B_value;
     C_tmp[1][1] += A_ptr[1] * B_value;
     C_tmp[2][1] += A_ptr[2] * B_value;
     C_tmp[3][1] += A_ptr[3] * B_value;

     B_value = B_ptr[2*ldb];
     C_tmp[0][2] += A_ptr[0] * B_value;
     C_tmp[1][2] += A_ptr[1] * B_value;
     C_tmp[2][2] += A_ptr[2] * B_value;
     C_tmp[3][2] += A_ptr[3] * B_value;

     B_value = B_ptr[3*ldb];
     C_tmp[0][3] += A_ptr[0] * B_value;
     C_tmp[1][3] += A_ptr[1] * B_value;
     C_tmp[2][3] += A_ptr[2] * B_value;
     C_tmp[3][3] += A_ptr[3] * B_value;

     A_ptr += lda;
     B_ptr ++;
    }

    for (int i = 0;i < M_UNROLL;i++)
    {
        for (int j = 0;j < N_UNROLL;j++)
        {
            C(i,j) += C_tmp[i][j];
        }
    }
}

void GemmUnrollReg( int m, int n, int k,
	      double *ap, int lda,
      	      double *bp, int ldb,
      	      double *cp, int ldc )
{
  int i, j;

  for ( j = 0;j < n;j += N_UNROLL)
  {
    for ( i = 0;i < m;i += M_UNROLL)
    {
      GemmUnrollRegKernel( k, &A( i,0 ), lda, &B( 0, j ), ldb, &C(i,j ), ldc );
    }
  }

  return;
}
  
