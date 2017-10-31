#define chi( i ) x[ (i)*incx ]   // map chi( i ) to array x 
#define psi( i ) x[ (i)*incy ]   // map psi( i ) to array y

void Dots( int n, double *x, int incx, double *y, int incy, double *gamma )
{
  int i;

  for ( i=0; i<n; i++ )
    *gamma = chi( i ) * psi( i ) + *gamma;
  
  return;
}
