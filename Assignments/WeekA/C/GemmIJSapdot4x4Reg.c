#define A( i,j ) *( ap + (j)*lda + i )   // map A( i,j ) to array ap in column-major order
#define B( i,j ) *( bp + (j)*ldb + i )   // map B( i,j ) to array bp in column-major order
#define C( i,j ) *( cp + (j)*ldc + i )   // map C( i,j ) to array bp in column-major order

void Sapdot4x4( int, double *, int, double *, int, double *, int );
  
void GemmIJSapdot4x4Reg( int m, int n, int k,
		      double *ap, int lda,
		      double *bp, int ldb,
		      double *cp, int ldc )
{
  int i, j;

  for ( i=0; i<m; i+=4 )
    for ( j=0; j<n; j+=4 ){
      Sapdot4x4( k, &A( i,0 ), lda, &B( 0, j ), ldb, &C(i,j ), ldc );
    }
  return;
}


void Sapdot4x4( int k,
		double *A, int lda,
		double *B, int ldb,
		double *C, int ldc )
{
  int p;
  register double
    c00=0.0, c01=0.0, c02=0.0, c03=0.0,
    c10=0.0, c11=0.0, c12=0.0, c13=0.0,
    c20=0.0, c21=0.0, c22=0.0, c23=0.0,
    c30=0.0, c31=0.0, c32=0.0, c33=0.0,
    areg,
    bp0, bp1, bp2, bp3;
  double *b0pntr=B, *b1pntr=B+ldb, *b2pntr=b1pntr+ldb, *b3pntr=b2pntr+ldb;

  for ( p=0; p<k; p++ ){
    areg = *A;
    
    bp0 = *b0pntr++;
    bp1 = *b1pntr++;
    bp2 = *b2pntr++;
    bp3 = *b3pntr++;
    
    c00 += areg * bp0; c01 += areg * bp1; c02 += areg * bp2; c03 += areg * bp3;
    areg = *(A+1);
    c10 += areg * bp0; c11 += areg * bp1; c12 += areg * bp2; c13 += areg * bp3;
    areg = *(A+2);
    c20 += areg * bp0; c21 += areg * bp1; c22 += areg * bp2; c23 += areg * bp3;
    areg = *(A+3);
    c30 += areg * bp0; c31 += areg * bp1; c32 += areg * bp2; c33 += areg * bp3; 

    A += lda;
}

  
  *C     += c00; *(C+ldc  ) += c01; *(C+2*ldc  ) += c02; *(C+3*ldc  ) += c03;
  *(C+1) += c10; *(C+ldc+1) += c11; *(C+2*ldc+1) += c12; *(C+3*ldc+1) += c13;
  *(C+2) += c20; *(C+ldc+2) += c21; *(C+2*ldc+2) += c22; *(C+3*ldc+2) += c23;
  *(C+3) += c30; *(C+ldc+3) += c31; *(C+2*ldc+3) += c32; *(C+3*ldc+3) += c33;    
}
    
