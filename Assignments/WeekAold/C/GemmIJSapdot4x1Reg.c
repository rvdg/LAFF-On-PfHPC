#define A( i,j ) *( ap + (j)*lda + i )   // map A( i,j ) to array ap in column-major order
#define B( i,j ) *( bp + (j)*ldb + i )   // map B( i,j ) to array bp in column-major order
#define C( i,j ) *( cp + (j)*ldc + i )   // map C( i,j ) to array bp in column-major order

void Sapdot4x1( int, double *, int, double *, double * );
  
void GemmIJSapdot4x1Reg( int m, int n, int k,
		      double *ap, int lda,
		      double *bp, int ldb,
		      double *cp, int ldc )
{
  int i, j;

  for ( i=0; i<m; i+=4 )
    for ( j=0; j<n; j++ ){
      Sapdot4x1( k, &A( i,0 ), lda, &B( 0, j ), &C(i,j ) );
    }
  return;
}


void Sapdot4x1( int n,
	     double *x, int incx,
	     double *y, 
	     double *gamma )
{
  int i;
  register double c0=0.0, c1=0.0, c2=0.0, c3=0.0, ycur;

  for ( i=0; i<n; i++ ){
    ycur = *y++;
    c0 += *x * ycur;
    c1 += *(x+1) * ycur;
    c2 += *(x+2) * ycur;
    c3 += *(x+3) * ycur;

    x += incx;
  }
  
  *gamma +=  c0;
  *(gamma+1) += c1;
  *(gamma+2) += c2;
  *(gamma+3) += c3;
}
    
